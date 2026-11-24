# AV3: Fourier series of $|x|^s$, $\operatorname{sgn}(x)|x|^s$ and their rotated mixtures, $s\in[2,4]$

**Bottom line.** The exponent homotopy does break parity, but $\zeta(3)$ enters this path only through the *odd* sector and only through *constant-sign linear functionals* of the coefficients: (i) the sum of sine coefficients $W(s)=\sum b_n(s)$ (the conjugate/Hilbert value at $0$, a log–sin integral), and (ii) the cross-parity pairing $B(s,s')=\sum a_n(s)b_n(s')$ (a log–cos integral). Every same-parity point value and Parseval identity of the family stays in the Bernoulli/Euler world; $\zeta(3),\zeta'(3),\zeta'(4),G,\beta'(3),\beta(4)$ appear in the *kink parts* of their $s$-derivatives and are exactly cancelled by the boundary tails (we exhibit each cancellation to 42–101 digits, and each cancellation is itself a proved closed-form corollary). The even family $|x|^s$ is proved to have **no** $1/n^3$ Fourier mode at *any* real $s$ — Aritro's Section-VI obstruction persists along the whole even path; in the odd family the kink amplitude $\propto\cos(\pi s/2)$ is *maximal* exactly at $s=2$, where the tail is $-\frac{4}{\pi}(1-(-1)^n)n^{-3}$.

---

## 1. Task 1 — Closed forms and asymptotics of the coefficients (PROVED, verified to 60 digits, step_05/06)

For $n\ge1$, real $s>-1$:
$$a_n(s)=\tfrac2\pi\!\int_0^\pi\! x^s\cos nx\,dx=\tfrac{2\pi^{s}}{s+1}\,{}_1F_2\Big(\tfrac{s+1}{2};\tfrac12,\tfrac{s+3}{2};-\tfrac{n^2\pi^2}{4}\Big)=\Re\big[\tfrac2\pi(-in)^{-s-1}\gamma(s+1,-in\pi)\big],$$
$$b_n(s)=\tfrac2\pi\!\int_0^\pi\! x^s\sin nx\,dx=\tfrac{2n\pi^{s+1}}{s+2}\,{}_1F_2\Big(\tfrac{s+2}{2};\tfrac32,\tfrac{s+4}{2};-\tfrac{n^2\pi^2}{4}\Big)=\Im\big[\tfrac2\pi(-in)^{-s-1}\gamma(s+1,-in\pi)\big],$$
$a_0(s)=2\pi^s/(s+1)$. The ${}_1F_2$ is a *generalized* Lommel function (its upper parameter $\ne1$); at $s=1$ it collapses to the elementary Lommel $s_{-1/2,3/2}(z)=(1-\cos z-z\sin z)\,z^{-3/2}$ (verified to 102 digits, step_48). All three forms agree with quadrature to 60 digits.

**Asymptotics** (from $\gamma=\Gamma-\Gamma_{\rm upper}$ + the large-$z$ expansion; verified in step_06, errors decay with boundary order $J$):
$$a_n(s)\sim \underbrace{-\tfrac2\pi\Gamma(s+1)\sin\tfrac{\pi s}{2}\,n^{-(s+1)}}_{\text{kink at }0}+\tfrac{2(-1)^n}{\pi}\sum_{j\ge1}(-1)^{j-1}\tfrac{\Gamma(s+1)}{\Gamma(s-2j+2)}\pi^{s-2j+1}n^{-2j},$$
$$b_n(s)\sim \underbrace{\tfrac2\pi\Gamma(s+1)\cos\tfrac{\pi s}{2}\,n^{-(s+1)}}_{\text{kink at }0}+\tfrac{2(-1)^{n+1}}{\pi}\sum_{j\ge0}(-1)^{j}\tfrac{\Gamma(s+1)}{\Gamma(s-2j+1)}\pi^{s-2j}n^{-(2j+1)}.$$
At integer $s$ the boundary series **terminate** ($1/\Gamma$ hits poles), giving exact finite forms (verified, step_06): e.g. $a_n(2)=4(-1)^n/n^2$, $a_n(3)=6\pi(-1)^n/n^2+12(1-(-1)^n)/(\pi n^4)$, $a_n(4)=(-1)^n(8\pi^2/n^2-48/n^4)$; $b_n(1)=2(-1)^{n+1}/n$, $b_n(2)=2\pi(-1)^{n+1}/n-4(1-(-1)^n)/(\pi n^3)$, $b_n(3)=2\pi^2(-1)^{n+1}/n+12(-1)^n/n^3$, $b_n(4)=2\pi^3(-1)^{n+1}/n+24\pi(-1)^n/n^3+48(1-(-1)^n)/(\pi n^5)$.

**The $1/n^3$ census (answer to Task 1's question).**
- **$f$-family (even): never.** The kink carries $1/n^3$ only when $s+1=3$, i.e. $s=2$, but its amplitude $\sin(\pi s/2)$ has a **node exactly there**; the boundary series runs over *even* powers $n^{-2j}$ only. So no $|x|^s$ has any $1/n^3$ mode at any real $s$.
- **$g$-family (odd): at $s=2$ with amplitude $-4/\pi$.** Net $1/n^3$ coefficient $=-\frac4\pi(1-(-1)^n)$: it lives entirely on **odd $n$**. Partner sums: weight $1\to\zeta(3)$; weight $(-1)^n\to-\eta(3)=-\frac34\zeta(3)$; odd-$n\to\frac78\zeta(3)$; weight $\sin(n\pi/2)\to\beta(3)=\pi^3/32$ (elementary — Euler!); root-of-unity weights $\to$ Hurwitz-$\zeta(3,p/q)$ combinations (not elementary). At $s=3$ the kink is off and the boundary gives $12(-1)^n/n^3$; at $s=4$, $24\pi(-1)^n/n^3$. Kink-amplitude census verified numerically (steps_38/40): e.g. $n^{3}b_n(2)\to-4/\pi=-1.273239545$ exactly as predicted.

---

## 2. Task 2 — The identity family, expansions, and where $\zeta(3)$ is

### 2.1 The eight natural identities (all PROVED exact; series side verified 13–102 digits)

| name | series | closed form | source |
|---|---|---|---|
| $Q(s)$ | $\sum a_n(s)$ (value at $0$) | $-\pi^s/(s+1)$ | Dirichlet |
| $P(s)$ | $\sum(-1)^n a_n(s)$ (at $\pi$) | $s\pi^s/(s+1)$ | Dirichlet |
| $R(s)$ | $\sum a_n(s)\cos(n\pi/2)$ (at $\pi/2$) | $(\pi/2)^s-\pi^s/(s+1)$ | Dirichlet |
| $\Phi(s)$ | $\sum a_n(s)^2$ (Parseval) | $2\pi^{2s}\big[\frac1{2s+1}-\frac1{(s+1)^2}\big]$ | Parseval |
| $T(s)$ | $\sum b_n(s)\sin(n\pi/2)$ | $(\pi/2)^s$ | Dirichlet |
| $\Psi(s)$ | $\sum b_n(s)^2$ | $2\pi^{2s}/(2s+1)$ | Parseval |
| $V(s)$ | $\sum b_n(s)/n$ | $\pi^{s+1}/((s+1)(s+2))$ | pair with $(\pi-x)/2$ |
| $W(s)$ | $\sum b_n(s)$ (conjugate at $0$) | $\frac1\pi\int_0^\pi t^s\cot(t/2)dt$ = log-sin values at integer $s$ | Hilbert |
plus mixed Parsevals $\sum a_n(s)a_n(s')=2\pi^{s+s'}\big[\frac1{s+s'+1}-\frac1{(s+1)(s'+1)}\big]$, $\sum b_n(s)b_n(s')=2\pi^{s+s'}/(s+s'+1)$.

Note $Q,P,R,\Phi,T,\Psi,V$ are **elementary for all $s$** (spot-checked at $s=2.5$, step_51). Corollaries at integer points: $T(2)$ reproves $\beta(3)=\pi^3/32$, $T(4)$ gives $\beta(5)=5\pi^5/1536$ — Euler's $\beta(2k+1)$ values come out of this family for free.

### 2.2 Expansions at $s=2,3,4$ (two $s$-derivatives, exact, from step_41/tables.json; series side spot-verified to 12–101 digits in steps_23/24/42/45)

With $L=\log\pi$, $\ell=\log(\pi/2)$:

| $F$ | $s$ | $F(s)$ | $F'(s)$ | $F''(s)$ |
|---|---|---|---|---|
| $Q$ | 2 | $-\pi^2/3$ | $\frac{\pi^2}{9}(1-3L)$ | $\frac{\pi^2}{27}(-9L^2+6L-2)$ |
| $Q$ | 3 | $-\pi^3/4$ | $\frac{\pi^3}{16}(1-4L)$ | $\frac{\pi^3}{32}(-8L^2+4L-1)$ |
| $Q$ | 4 | $-\pi^4/5$ | $\frac{\pi^4}{25}(1-5L)$ | $\frac{\pi^4}{125}(-25L^2+10L-2)$ |
| $P$ | 2 | $2\pi^2/3$ | $\frac{\pi^2}{9}(6L+1)$ | $\frac{2\pi^2}{27}(9L^2+3L-1)$ |
| $P$ | 3 | $3\pi^3/4$ | $\frac{\pi^3}{16}(12L+1)$ | $\frac{\pi^3}{32}(24L^2+4L-1)$ |
| $P$ | 4 | $4\pi^4/5$ | $\frac{\pi^4}{25}(20L+1)$ | $\frac{2\pi^4}{125}(50L^2+5L-1)$ |
| $R$ | 2 | $-\pi^2/12$ | $\frac{\pi^2}{36}(4-3L-9\log2)$ | $\frac{\pi^2}{4}\ell^2-\frac{\pi^2}{3}L^2+\frac{2\pi^2}{9}L-\frac{2\pi^2}{27}$ |
| $R$ | 3 | $-\pi^3/8$ | $\frac{\pi^3}{16}(1-2L-2\log2)$ | $\frac{\pi^3}{32}(4\ell^2-8L^2+4L-1)$ |
| $R$ | 4 | $-11\pi^4/80$ | $\frac{\pi^4}{400}(16-55L-80\log2)$ | $\frac{\pi^4}{16}\ell^2-\frac{\pi^4}{5}L^2+\frac{2\pi^4}{25}L-\frac{2\pi^4}{125}$ |
| $\Phi$ | 2 | $8\pi^4/45$ | $\frac{8\pi^4}{675}(30L-1)$ | $\frac{4\pi^4}{3375}(600L^2-40L-17)$ |
| $\Phi$ | 3 | $9\pi^6/56$ | $\frac{3\pi^6}{784}(84L-5)$ | $\frac{\pi^6}{21952}(14112L^2-1680L-5)$ |
| $\Phi$ | 4 | $32\pi^8/225$ | $\frac{16\pi^8}{10125}(180L-11)$ | $\frac{4\pi^8}{455625}(64800L^2-7920L+313)$ |
| $T$ | 2 | $\pi^2/4$ | $\frac{\pi^2}{4}\ell$ | $\frac{\pi^2}{4}\ell^2$ |
| $T$ | 3 | $\pi^3/8$ | $\frac{\pi^3}{8}\ell$ | $\frac{\pi^3}{8}\ell^2$ |
| $T$ | 4 | $\pi^4/16$ | $\frac{\pi^4}{16}\ell$ | $\frac{\pi^4}{16}\ell^2$ |
| $\Psi$ | 2 | $2\pi^4/5$ | $\frac{4\pi^4}{25}(5L-1)$ | $\frac{8\pi^4}{125}(25L^2-10L+2)$ |
| $\Psi$ | 3 | $2\pi^6/7$ | $\frac{4\pi^6}{49}(7L-1)$ | $\frac{8\pi^6}{343}(49L^2-14L+2)$ |
| $\Psi$ | 4 | $2\pi^8/9$ | $\frac{4\pi^8}{81}(9L-1)$ | $\frac{8\pi^8}{729}(81L^2-18L+2)$ |
| $V$ | 2 | $\pi^3/12$ | $\frac{\pi^3}{144}(12L-7)$ | $\frac{\pi^3}{1728}(144L^2-168L+74)$ |
| $V$ | 3 | $\pi^4/20$ | $\frac{\pi^4}{400}(20L-9)$ | $\frac{\pi^4}{8000}(400L^2-360L+122)$ |
| $V$ | 4 | $\pi^5/30$ | $\frac{\pi^5}{900}(30L-11)$ | $\frac{\pi^5}{27000}(900L^2-660L+182)$ |

All elementary — **no odd $\zeta$ in the value or in any $s$-derivative of these seven identities, at any $s$**.

### 2.3 The kink-split analysis: $\zeta(3),\zeta'(3),\zeta'(4),G,\beta'(3),\beta(4)$ all appear — and cancel (each line verified to the stated digits)

Splitting the series side of an identity into kink + boundary parts, the odd-weight constants enter via the kink:
- $Q'(2)$: kink part $=\frac{d}{ds}\big[-\frac2\pi\Gamma(s+1)\sin\frac{\pi s}{2}\zeta(s+1)\big]_{s=2}=2\zeta(3)$; verified $\sum_n\big[a_n'(2)-2/n^3\big]=Q'(2)-2\zeta(3)$ to **101 digits** (step_23).
- $Q''(2)$: kink part $=4\zeta'(3)$; verified $\sum_n\big[a_n''(2)+4\log n/n^3\big]=Q''(2)-4\zeta'(3)$ to **42 digits** (step_24).
- $V'(1)$: kink part $=-\zeta(3)$ (unit coefficient); verified to **101 digits** (step_23).
- $T'(1)$: kink part $=-G$ (Catalan); verified to **~100 digits** (step_48).
- $T'(2)$: kink part $=-\frac4\pi(\psi(3)\beta(3)+\beta'(3))$; verified to **101 digits** (step_29).
- $T'(3)$: kink part $=6\beta(4)$; verified, residual $0.0$ (step_38).
- $W'(2)$: kink part $=-\frac4\pi(\zeta'(3)+\psi(3)\zeta(3))$; verified, residual $0.0$ (step_48).

Each cancellation is a theorem: the exact closed form forces the boundary tail to eat the kink constant. The obstruction to $\zeta(3)$ is *precisely* the rationality of $\zeta(2n)/\pi^{2n}$.

### 2.4 Corollary constants (all PROVED via the identities named; all verified)

| constant | value | proved via | verified |
|---|---|---|---|
| $\sum \operatorname{Si}(n\pi)/n^3$ | $5\pi^3/72$ | $Q'(2)$ | 90 digits, 3 methods |
| $\sum \operatorname{Si}_2(n\pi)/n^3$, $\operatorname{Si}_2(z)=\int_0^z\frac{\operatorname{Si}u}{u}du$ | $49\pi^3/432$ | $Q''(2)$ | 15 digits + proof |
| $\sum \operatorname{Si}(n\pi)\sin(n\pi/2)/n^2$ | $\pi^2\log2/4$ | $T'(1)$ | full precision |
| $\sum(\operatorname{Si}(n\pi)-\frac\pi2)\sin(n\pi/2)/n^4$ | $\frac{\pi^4(5-2\log2)}{192}-\frac\pi2\beta(4)$ | $T'(3)$ | 101 digits |
| $\sum_k(-1)^k\operatorname{Ci}((2k+1)\pi)/(2k+1)^3$ | $\frac{\pi^3}{32}(\log(4\pi)-3+\gamma)-\beta'(3)$ | $T'(2)$ | 101 digits |
| $\sum_n(-1)^n\operatorname{Ci}(n\pi)/n^4$ | $\frac{\pi^4}{2880}(45+4\log2-28\gamma-28\log\pi)+\frac78\zeta'(4)$ | $P'(3)$ | 101 digits |
| $\operatorname{Si}_2(z)-\frac\pi2\log z\to C$ | $C=\frac{\pi\gamma}{2}$ | $\Gamma(s)\sin\frac{\pi s}{2}\big|_{s=0}'$ | 40 digits |

### 2.5 The identities that DO produce $\zeta(3)$ (PROVED, verified 3 ways to 80 digits)

**$W(s)=\sum_n b_n(s)=\frac1\pi\int_0^\pi t^s\cot(t/2)dt$** (sum of sine coefficients = conjugate value at 0):
$$W(1)=2\log2,\qquad \boxed{W(2)=2\pi\log2-\frac{7\zeta(3)}{\pi}},\qquad W(3)=2\pi^2\log2-9\zeta(3),\qquad W(4)=2\pi^3\log2-18\pi\zeta(3)+\frac{93\zeta(5)}{\pi}.$$
Its kink part is $W_K(s)=\frac2\pi\Gamma(s+1)\cos\frac{\pi s}{2}\zeta(s+1)$: at $s=2$ this is $-\frac4\pi\zeta(3)$ with the antinode (maximal amplitude), and it **survives** because the functional has constant sign on the tail. This is the family member that plays the role of "Parseval for $x$ gives $\zeta(2)$". The full odd tower $\zeta(3),\zeta(5),\dots,\zeta(2k+1)$ appears at $s=2k$.

**Cross-parity pairing $B(s,s')=\sum_n a_n(s)b_n(s')=-\frac1\pi\langle f_s,\tilde g_{s'}\rangle$** (log-cos integrals; e.g. $B(1,1)=\frac4\pi\int_0^\pi x\log(2\cos(x/2))dx$, verified to 36 digits):
$$B(1,1)=-\frac{7\zeta(3)}{\pi},\quad B(2,1)=-8\zeta(3),\quad B(1,2)=-7\zeta(3)+\frac{31\zeta(5)}{\pi^2},\quad B(2,2)=-8\pi\zeta(3)+\frac{31\zeta(5)}{\pi},$$
$$B(1,3)=-7\pi\zeta(3)+\frac{93\zeta(5)}{2\pi},\ B(2,3)=-8\pi^2\zeta(3)+48\zeta(5),\ B(1,4)=-7\pi^2\zeta(3)+93\zeta(5)-\frac{381\zeta(7)}{\pi^2},\ B(2,4)=-8\pi^3\zeta(3)+96\pi\zeta(5)-\frac{381\zeta(7)}{\pi},$$
$$B(3,3)=-12\pi^3\zeta(3)+\frac{237\pi\zeta(5)}2-\frac{1143\zeta(7)}{4\pi},\ B(3,4)=-12\pi^4\zeta(3)+\frac{381\pi^2\zeta(5)}2-1143\zeta(7)+\frac{4599\zeta(9)}{2\pi^2},\ B(4,4)=-16\pi^5\zeta(3)+288\pi^3\zeta(5)-1914\pi\zeta(7)+\frac{4599\zeta(9)}{\pi}.$$

### 2.6 OPEN constants (high precision; PSLQ null at two+ precisions)

- $\mathrm{CiS}=\sum_n\operatorname{Ci}(n\pi)/n^3=0.07116843744486633624888448010900074229067087805042$ (103 digits, two methods agreeing to $8.9\times10^{-103}$). PSLQ against 15 named constants at 100–110 digits: only precision-floor relations (residual $\sim10^{-99}$, maxcoeff $\sim2\times10^9$; residual tracks the floor at 50/70/90/110 dps). **OPEN.**
- $W'(2)=2\pi\log\pi\log2-\frac4\pi\big(\mathrm{CiS}+\zeta'(3)+\tfrac34\log\pi\,\zeta(3)-\gamma\zeta(3)\big)-\frac{21\zeta(3)}{2\pi}=0.69898415458911186737886913136066417279823793041123$ (three methods, 45 digits; the reduction to CiS is PROVED). **OPEN.**
- $R_4=\sum_n(\operatorname{Si}(n\pi)-\frac\pi2)/n^4=0.27264972611053263981893507622668795678657010150647$; **OPEN.**
- $W'(3)=2\pi^2\log\pi\log2-(9\log\pi+\tfrac{15}{2})\zeta(3)+\frac{\pi^4}{15}+\frac{12}{\pi}R_4=1.79808054679361434950988293086017951541326397$ (three methods, 45 digits; reduction PROVED). **OPEN.**

---

## 3. Task 3 — The rotated family $h_{s,\theta}=\cos\theta\,f_s+\sin\theta\,g_s$ (verified, step_47/48)

- **Parseval is exactly $\theta$-blind**: $\frac1\pi\int_{-\pi}^{\pi}h_{s,\theta}^2=\cos^2\theta\,\Phi_{\rm full}+\sin^2\theta\,\Psi=\frac{2\pi^{2s}}{2s+1}$ for every $\theta$ ($\int f_s^2=\int g_s^2$; cross terms die by parity). $L^2$ never sees the rotation.
- **Point values are tautologies for every $\theta$** (series = function value, verified to $2\times10^{-81}$); at integer $s$ every point evaluation splits into elementary (sawtooth/Bernoulli) pieces only.
- **The $\zeta(3)$ content scales with $\sin\theta$** (conjugate-sum/$W$-type): $\sum_n d_n(\theta)=\sin\theta\big(2\pi\log2-\frac{7\zeta(3)}{\pi}\big)$ at $s=2$ — extracting $\zeta(3)$ exactly at $\theta=\pi/6,\pi/4,\pi/3,\pi/2$ (verified). Cleanest at **$\theta=\pi/2$** (pure odd = maximal parity break).
- **Cross-pairing** $\sum c_nd_n=\sin\theta\cos\theta\,B(s,s)=\tfrac12\sin(2\theta)\,B(s,s)$: maximal at **$\theta=\pi/4$** (e.g. $-\frac{7\zeta(3)}{2\pi}$ at $s=1,\theta=\pi/4$; partial sums to $8\times10^{-12}$).
- So: *no point makes the parity break produce $\zeta(3)$*; the parity-breaking **linear** functionals do — at $\theta=\pi/2$ (W) or $\theta=\pi/4$ (B).

## 4. Task 4 — Plain words, and the squircle comparison

The linear homotopy $(1-t)x^2+tx^4$ stays even, and we proved the even family has *no* $1/n^3$ mode at any real exponent: the kink amplitude $\sin(\pi s/2)$ has its zero exactly where $n^{-(s+1)}=n^{-3}$. The non-linear exponent path *does* rotate through the $1/n^3$ channel, but in the even sector it passes through a **node**; only in the odd sector ($g_s$, or $\theta=\pi/2$ in the rotated family) does the $1/n^3$ tail have maximal, nonzero amplitude — and it lives on odd $n$, so it is invisible to point evaluations (alternating weights, all tautologies anyway) and to Parseval (quadratic: pairs $n^{-a}n^{-b}$ into even or $\eta$-chains, i.e. $\zeta(2k),\zeta'(2k)$ only). $\zeta(3)$ survives precisely under **constant-sign summation of the tail**: the sum of all sine coefficients (= the conjugate function at 0 = a log-sin/Clausen integral) and the cos×sin pairing (= a log-cos integral). In the $s$-derivatives of the elementary identities, $\zeta(3),\zeta'(3),\zeta'(4),G,\beta'(3),\beta(4)$ are all present "virtually" in the kink parts and are exactly cancelled by the boundary tails — each cancellation being one of our proved Si/Ci corollaries. Versus the squircle: there $\zeta(3)$ also appears at an **endpoint** of the shape path ($p=1$, diamond), also via a boundary-layer object (the corner coefficients $M_k$), also partnered with $\log2$; here the endpoint is $s=2$ in the odd sector, the boundary layer is the kink at $0$, and the certificate is again a Clausen-type integral. Difference: the squircle's corner constants are positive finite combinations, while this path yields $\zeta(3)$ only through signed infinite sums that need the kink/boundary split — and (caveat for the paper's program) the Fourier route repackages the classical log-sin proof of $\zeta(3)$ rather than avoiding it.

## 5. Evidence table (constant — value — status — digits/methods — step files)

| quantity | value | status | evidence | steps |
|---|---|---|---|---|
| coeff closed forms (${}_1F_2$, $\gamma$) | §1 | PROVED | 60 digits vs quadrature | 05 |
| kink+boundary asymptotics | §1 | PROVED | errors decay with $J$, table | 06, 40 |
| exact integer forms $s=1..4$ | §1 | PROVED (IBP) | consistent everywhere | 06 |
| $W(2)=2\pi\log2-7\zeta(3)/\pi$ | 1.676786008389011865250305814789164935495 | PROVED | 80 digits, 3 methods | 22 |
| $W(3),W(4)$ | 2.863664799277884521…, 5.705095152293848340… | PROVED | 80 digits, 3 methods | 22 |
| $B(s,s')$ tower (11 values) | §2.5 | PROVED | partial sums $10^{-7}$–$10^{-12}$; $B(1,1)$ integral 36 digits | 34–37 |
| kink splits $Q'(2),V'(1),T'(1),T'(2),T'(3),W'(2)$ | §2.3 | PROVED | 100–101 digits ($Q''(2)$: 42) | 23–29, 38, 48 |
| $\sum\mathrm{Si}(n\pi)/n^3=5\pi^3/72$ | 2.1532136583541541788525218796598191112656450392976 | PROVED via $Q'(2)$ | 90 digits, 3 methods | 08, 14 |
| $\sum\mathrm{Si}_2(n\pi)/n^3=49\pi^3/432$ | 3.5169156419784518254591190701110378817338869 | PROVED via $Q''(2)$ | 15 digits + proof | 12–19 |
| $\sum\mathrm{Si}(n\pi)\sin(n\pi/2)/n^2=\pi^2\log2/4$ | 1.71027211596427913621186978848852401782494 | PROVED via $T'(1)$ | full precision | 26 |
| $\mathrm{S}_4=\sum(\mathrm{Si}-\frac\pi2)\sin(n\pi/2)/n^4$ | $\frac{\pi^4(5-2\log2)}{192}-\frac\pi2\beta(4)=0.279943392949165487288669609612754993174355$ | PROVED via $T'(3)$ | 101 digits | 38, 39 |
| $\mathrm{Cisum}=\sum_k(-1)^k\mathrm{Ci}((2k+1)\pi)/(2k+1)^3$ | $\frac{\pi^3}{32}(\log(4\pi)-3+\gamma)-\beta'(3)=0.0733015660215609599845396324296814531596501$ | PROVED via $T'(2)$ | 101 digits | 26, 29 |
| $\mathrm{CiS4}=\sum(-1)^n\mathrm{Ci}(n\pi)/n^4$ | $\frac{\pi^4}{2880}(45+4\log2-28\gamma-28\log\pi)+\frac78\zeta'(4)=-0.0752431151829964334833783188165789047850027$ | PROVED via $P'(3)$ | 101 digits | 43–45 |
| $\mathrm{CiS}=\sum\mathrm{Ci}(n\pi)/n^3$ | 0.07116843744486633624888448010900074229067087805042 | **OPEN** (PSLQ null, 4 precisions) | 103 digits, 2 methods | 30–33 |
| $W'(2)$ | 0.69898415458911186737886913136066417279823793041123 | reduction PROVED; CiS OPEN | 45 digits, 3 methods | 30–33, 48 |
| $R_4,\ W'(3)$ | 0.27264972611053263981893507622668795678657010150647; 1.79808054679361434950988293086017951541326397 | reduction PROVED; $R_4$ OPEN | 45 digits, 3 methods; PSLQ spurious (floor test) | 49, 50 |
| rotation laws | $\sin\theta\cdot W$, $\tfrac12\sin2\theta\cdot B$, Parseval $\theta$-blind | PROVED | $\zeta(3)$ extracted exactly at 4 angles; point value $2\times10^{-81}$ | 47, 48 |
| expansion tables $F,F',F''$ at $s=2,3,4$ | §2.2 (tables.json) | PROVED | series-side 12–101 digits | 41–45 |

**Files** (in `new_results/AV3_fractional_fourier`): `step_01.py`–`step_59.py` (every computation; numbering follows the run order — key ones: `step_05/06` coefficients & asymptotics, `step_07` derivative formulas, `step_08/14` Si-identity, `step_11` $C=\pi\gamma/2$, `step_19` $\mathrm{Si}_2$, `step_22` W-family, `step_23–29` kink splits & T-identities, `step_33/50` PSLQ-null tests, `step_35/37` B-tower, `step_38/39` $T'(3)$ & $\mathrm{S}_4$, `step_41` tables, `step_42–45` s=3,4 verifications & CiS4, `step_47/48` rotation & Lommel, `step_49` $W'(3)/R_4$), `progress.log` (one line per finished step), `tables.json` (exact expansion tables), `B_forms.json` (B-tower closed forms), `humps.json`, `h_coeffs.json` (Si₂ asymptotics data).