# D12 — Asymptotic/Taylor expansions of the Adler–Tanton / Keller–Vakil constant $\pi_p$

All numbers below were printed by the code listed at the end (run in this session at 60–150 digit precision). Labels: **PROVED** / **NUMERICAL** / **CONJECTURE**.

---

## 1. Statements

### 1.1 The integral and two exact representations (PROVED)

From $\pi_p=\frac2p\int_0^1\!\big(t^{1-p}+(1-t)^{1-p}\big)^{1/p}dt$, symmetry $t\leftrightarrow1-t$ plus factoring out $t^{-(p-1)}$ on $[0,\tfrac12]$ gives, with $h=1/p$,

$$\pi_p=\frac{4}{p}\int_0^{1/2}t^{-1+1/p}\Big(1+\big(\tfrac{t}{1-t}\big)^{p-1}\Big)^{1/p}dt
=2h\int_0^\infty\cosh\!\Big(\frac u2\Big)^{-1-h}\cosh\!\Big(\frac{u(1-h)}{2h}\Big)^{h}du .$$

**Exact splitting (the backbone of everything):**
$$\boxed{\;\pi_p=4\cdot 2^{-1/p}+R(p)\;}\qquad
R(p)=\frac{4}{p(p-1)}\int_0^\infty e^{-\frac{s}{p(p-1)}}\Big(1+e^{-\frac{s}{p-1}}\Big)^{-1-\frac1p}\Big[\big(1+e^{-s}\big)^{1/p}-1\Big]ds .$$

(The "$4\cdot2^{-1/p}$" piece is also $g(h)=4h\int_0^1v^{h-1}(1+v)^{-1-h}dv=4\,{}_2F_1(1+h,h;1+h;-1)=4\cdot2^{-h}$, since ${}_2F_1(a,b;a;z)=(1-z)^{-b}$.)

### 1.2 Formula and duality verified (NUMERICAL, 40–80 digits)

- $\pi_2=\pi$ (|diff| $2.2\times10^{-37}$), $\pi_1=4$ exactly, $\pi_{10^6}=3.99999722741224\to4=\pi_\infty$.
- Duality $\pi_p=\pi_q$, $1/p+1/q=1$: pairs $(1.25,5),(1.5,3),(2.5,\tfrac53),(3,1.5),(4,\tfrac43),(6,1.2),(10,\tfrac{10}9)$ agree to $\le 3.6\times10^{-71}$ (exact conjugates) — see evidence table.
- Minimum at $p=2$ confirmed (values at $1.9,1.95,2,2.05,2.1$).

### 1.3 Expansion at $p=\infty$ (PROVED for $N_0..N_3$; the rest verified $\ge$100 digits)

$\pi_p=4\big(1+\frac{N_0}{p}+\frac{N_1}{p^2}+\frac{N_2}{p^3}+\frac{N_3}{p^4}+\frac{N_4}{p^5}+\frac{N_5}{p^6}+O(p^{-7})\big)$ with

$$N_0=-\log2,\qquad N_1=\tfrac12\log^2 2,$$
$$N_2=\frac{\pi^2}{24}-\frac{\log^3 2}{6},\qquad
N_3=\frac{\log^4 2}{24}+\frac{\pi^2}{24}(1-\log2)+\frac{\zeta(3)}{4},$$
$$N_4=-\frac{\log^5 2}{120}+\frac{\pi^2}{24}\Big(1-\log2+\frac{\log^2 2}{2}\Big)+\frac{\zeta(3)}{4}(1-\log2)+\frac{\pi^4}{2880},$$
$$N_5=\frac{\log^6 2}{720}+\frac{\pi^2}{24}\Big(1-\log2+\frac{\log^22}{2}-\frac{\log^32}{6}\Big)+\frac{\zeta(3)}{4}\Big(1-\log2+\frac{\log^22}{2}\Big)+\frac{\pi^4}{2880}(1-\log2)-\frac{7\pi^4}{1920}-\frac{\zeta(3)\pi^2}{48}+\frac{\zeta(5)}{8}.$$

Equivalently $R(p)=\frac{\pi^2}{6}p^{-3}+\big(\frac{\pi^2}{6}+\zeta(3)-\frac{\pi^2\log2}{6}\big)p^{-4}+O(p^{-5})$ (**PROVED**, §2).
PSLQ at 95 digits: $24N_2=\pi^2-4\log^32$; $24N_3=\log^42+\pi^2-\pi^2\log2+6\zeta(3)$; $2880N_4=-24\log^52+120\pi^2(1-\log2+\frac{\log^22}2)+720\zeta(3)(1-\log2)+\pi^4$ (residuals $\sim10^{-97}$). **No $\sqrt2$, no $\log(1+\sqrt2)$, no $\operatorname{Li}_n(\sqrt2-1)$** (joint-basis PSLQ puts coefficient 0 on all Euclidean constants). Constants are $\{\log2,\pi^2,\zeta(3),\zeta(5)\}$ here; in general alternating Euler sums / $\operatorname{Li}_n(1/2)$ (see $J$-values in §2).

**Bonus structure (NUMERICAL, 100 digits):** writing $\pi_{1/h}=4\cdot2^{-h}\big(1+\sum_{m\ge3}\frac{a_m}{4}h^m\big)$,
$a_3=\frac{\pi^2}{6},\ a_4=\frac{\pi^2}{6}+\zeta(3),\ a_5=\frac{\pi^2}{6}+\zeta(3)+\frac{\pi^4}{720},\ a_6=\frac{\pi^2}{6}+\zeta(3)+\frac{\pi^4}{720}-\frac{7\pi^4}{480}-\frac{\zeta(3)\pi^2}{12}+\frac{\zeta(5)}{2}$.

### 1.4 Radius of convergence (NUMERICAL, with PROVED exact coefficient formula)

**Exact formula (PROVED):** $N_k=\dfrac{(-\log2)^{k+1}}{(k+1)!}+\displaystyle\sum_{m=0}^{k-1}\int_0^\infty[h^m]\Psi(s)\,ds$, where $[h^m]\Psi$ is an explicit **rational** polynomial in $s$, $L=\log(1+e^{-s})$, $\log2$; hence every $N_k$ is an exact rational combination of $\log2$ and $J_{j,i}=\int_0^\infty s^i\log^j(1+e^{-s})ds$ (alternating MZVs of weight $i+j+1$).

**Result: radius $=0$ — the series is divergent asymptotic.** $|N_k|$ grows factorially, $|N_k|\asymp \frac{2(k-2)!}{\pi^{k-1}}$ with a period-4 sign pattern:
- $|N_k|^{1/k}\approx k/(e\pi)$ (at $k=40$: $4.24$ vs $4.68$);
- $\big|N_{k+2}/N_k\big|^{1/2}\cdot\frac{\pi}{k+1}\to1$ (0.746 at $k=10$ → 0.943 at $k=40$);
- odd $k$: $N_k\pi^{k-1}/\big(2(k-2)!\big)\to(-1)^{(k-1)/2}$ (1.283 at $k=21$ → 1.124 at $k=39$, slow $1+O(1/k)$ approach).

Mechanism (heuristic, validated numerically): the integrand's branch points $h^*(s)=\frac{-i\pi(2\nu+1)}{s-i\pi(2\nu+1)}$ accumulate at $h=0$ along $\arg h=\pm\frac\pi2$; Cauchy integration confirms $[h^m]\psi_2(s)\approx-\frac2s\operatorname{Re}(1+is/\pi)^{m+1}$ (ratios 0.80–0.93 → 1). Hence Borel singularities at $\pm i\pi$, nonperturbative scale $e^{-\pi/h}$: optimal truncation at $p=10$ stops at $k=31\approx\pi/h=31.4$ with error $7.9\times10^{-17}\sim e^{-\pi/h}$.

### 1.5 Expansion at $p=1$ (via duality; NUMERICAL)

With the task's variable $\eta=\frac1q-1$ ($q\to1$): since $h=1/p=-\eta$,
$$\pi_q=4\Big(1+\sum_{k\ge1}(-1)^k N_{k-1}\,\eta^k\Big),\qquad\text{radius }0.$$
In the alternative variable $\eta=q-1$ ($h=\eta/(1+\eta)$): $\pi_{1+\eta}=4\sum M_n\eta^n$, $M_0=1$, $M_n=\sum_{j=1}^n(-1)^{n-j}\binom{n-1}{j-1}N_{j-1}$; e.g. $M_1=-\log2$, $M_2=\log2+\frac{\log^22}2$, $M_3=-\log2-\log^22+\frac{\pi^2}{24}-\frac{\log^32}6$; $|M_n|^{1/n}$ grows (3.75 at $n=39$ vs $n/(e\pi)=4.2$) ⇒ **radius 0 here too**. Verified against direct evaluation: at $\eta=0.02,0.05,0.1$ best-truncation errors $1.8\times10^{-45},8.8\times10^{-29},1.1\times10^{-16}$ at $n=40,40,34$ ($\pi/h=160,66,34.6$ ✓).

### 1.6 Comparison with the Euclidean case (Task 4)

- **Shared constants:** $\log 2$ (weight 1; our $N_0=-\log2$ is pure, Euclidean $M_0=\log2-\sqrt2\log(1+\sqrt2)$), $\pi^2$ (our $N_2,N_3,\dots$; Euclidean via the $\operatorname{Li}_2(\sqrt2-1)$ sector and the $\frac{\pi^2}{6}(p-1)^2$ term at the rhombus end), $\zeta(3)$ (our $N_3$, coefficient $+1/4$; Euclidean weight-3 $M_2$ sector).
- **Disjoint alphabets:** Euclidean lives in the $\mathbb Q(\sqrt2)$-alphabet $\{\sqrt2,\log(1+\sqrt2)=\operatorname{arsinh}1,\operatorname{Li}_n(\sqrt2-1)\}$; ours in the $\{\pm1\}$-colored (alternating) MZV alphabet $\{\log2,\eta(n),\operatorname{Li}_n(1/2),$ alternating double zetas$\}$. PSLQ offered both alphabets simultaneously puts weight only on ours.
- **Structural contrast:** the Euclidean $1/p$-series at $\infty$ **converges** ($p>1$) while its $p\to1$ series is conjectured to diverge. Here duality $\pi_p=\pi_q$ identifies the two ends by the Möbius map $h=\eta/(1+\eta)$, so **both** ends have radius 0 — one phenomenon, not two.

---

## 2. Proofs

**Lemma 1 (exact split).** Integrand symmetric under $t\mapsto1-t$, so $\pi_p=\frac4p\int_0^{1/2}$. On $[0,\frac12]$, $\big(t^{-(p-1)}+(1-t)^{-(p-1)}\big)^{1/p}=t^{-1+1/p}\big(1+(\frac{t}{1-t})^{p-1}\big)^{1/p}$. Splitting $\big(1+(\frac{t}{1-t})^{p-1}\big)^{1/p}=1+[\cdots]$ and using $\frac4p\int_0^{1/2}t^{-1+1/p}dt=4\cdot2^{-1/p}$ gives the split. $\square$

**Lemma 2 (change of variables).** $r=t/(1-t)$ sends $[0,\frac12]\to[0,1]$; then $r=e^{-s/(p-1)}$ gives the stated $R(p)$ (algebra: $r^{-1+1/p}dr=-\frac1{p-1}e^{-s/(p(p-1))}ds$). $\square$

**Lemma 3 (bounds; $h=1/p\le1/2$, $L=\log(1+e^{-s})$).** (i) $hL\le(1+e^{-s})^h-1\le hL+h^2L^2$ (from $x\le e^x-1\le x+x^2$ for $0\le x\le\log2$). (ii) $w(s,h):=e^{-sh^2/(1-h)}(1+e^{-sh/(1-h)})^{-1-h}$ satisfies $\tfrac14e^{-2sh^2}\le w\le1$. (iii) Uniformly for $0\le s\le h^{-1/2}$: $w(s,h)=\frac12\big[1+h\big(\frac s2-\log2\big)\big]+O\big((1+s^2)h^2\big)$. *Proof of (iii):* with $u=sh/(1-h)\le2\sqrt h$, $\log(1+e^{-u})=\log2-\frac u2+\frac{u^2}8+O(u^4)$ (converges for $|u|<\pi$), so $\log w=-sh^2/(1-h)-(1+h)\log(1+e^{-u})=-\log2+h(\frac s2-\log2)+O((1+s^2)h^2)$; exponentiate. $\square$

**Proposition.** $R(p)=\frac{\pi^2}{6}p^{-3}+\big(\frac{\pi^2}{6}+\zeta(3)-\frac{\pi^2\log2}{6}\big)p^{-4}+O(p^{-5})$.
*Proof.* By (i), $R=\frac{4h^2}{1-h}\big[hI_1+\frac{h^2}2I_2+O(h^3)\big]$, $I_m=\int_0^\infty wL^m\,ds$ (using $\int L^3<\infty$). By (iii), splitting at $s=h^{-1/2}$ with exponential tails: $I_1=\frac12J_{1,0}+\frac h2\big(\frac12J_{1,1}-\log2\,J_{1,0}\big)+O(h^2)$, $I_2=\frac12J_{2,0}+O(h)$. Since $\frac{4h^2}{1-h}=4h^2+4h^3+O(h^4)$: $R=2J_{1,0}h^3+4\big[\frac{J_{1,0}}2+\frac{J_{1,1}}4-\frac{\log2}2J_{1,0}+\frac{J_{2,0}}4\big]h^4+O(h^5)$. Substituting the evaluations below gives the claim. $\square$

**$J$-evaluations (PROVED).** $J_{1,k}=\int_0^\infty s^kL\,ds=\sum_{n\ge1}\frac{(-1)^{n-1}}n\int_0^\infty s^ke^{-ns}ds=k!\,\eta(k+2)$, so $J_{1,0}=\frac{\pi^2}{12}$, $J_{1,1}=\frac{3\zeta(3)}4$, $J_{1,2}=\frac{7\pi^4}{360}$, $J_{1,3}=\frac{45\zeta(5)}8$. For $J_{2,0}$: substitute $v=e^{-s}$, use $\log^2(1+v)=2\sum_{n\ge2}(-1)^n\frac{H_{n-1}}nv^n$ and the Euler sum $\sum_{n\ge1}(-1)^{n+1}\frac{H_n}{n^2}=\frac58\zeta(3)$: $J_{2,0}=2\sum(-1)^n\frac{H_{n-1}}{n^2}=\frac{\zeta(3)}4$. $\square$

**Corollary.** The formulas for $N_0,N_1,N_2,N_3$ follow from $\pi_p/4=2^{-h}+R/4$ and $2^{-h}=\sum\frac{(-\log2)^k}{k!}h^k$. $\square$

**$N_4,N_5$:** same machinery to order $h^5,h^6$ (coefficients $[h^m]\Psi$ computed exactly, §1.4), using the standard Euler-sum evaluations (each verified to 100 digits, provable by the same expansion/stuffle techniques):
$J_{3,0}=-\frac{\log^42}4+\frac{\pi^2\log^22}4+\frac{\pi^4}{15}-\frac{21}4\zeta(3)\log2-6\operatorname{Li}_4(\tfrac12)$,
$J_{2,1}=\frac{\log^42}6-\frac{\pi^2\log^22}6-\frac{\pi^4}{24}+\frac72\zeta(3)\log2+4\operatorname{Li}_4(\tfrac12)$,
$J_{4,0}=-\frac45\log^52+\frac23\pi^2\log^32-\frac{21}2\zeta(3)\log^22+24\zeta(5)-24\operatorname{Li}_4(\tfrac12)\log2-24\operatorname{Li}_5(\tfrac12)$,
$J_{3,1}=\frac25\log^52-\frac13\pi^2\log^32+\frac{21}4\zeta(3)\log^22-\frac12\zeta(3)\pi^2-\frac{99}{16}\zeta(5)+12\operatorname{Li}_4(\tfrac12)\log2+12\operatorname{Li}_5(\tfrac12)$,
$J_{2,2}=\frac13\zeta(3)\pi^2-\frac{29}8\zeta(5)$. Notable cancellations: $\frac1{12}J_{3,0}+\frac18J_{2,1}=\frac{\pi^4}{2880}$ and $\frac1{48}(J_{4,0}+2J_{3,1})=-\frac{\zeta(3)\pi^2}{48}+\frac{31\zeta(5)}{128}$. (Label: **PROVED modulo these standard Euler sums**, all verified $\ge100$ digits.)

---

## 3. Evidence table

**Duality (60 digits, cosh form):** $\pi_{1.25}=\pi_5=3.49845680804544907591825728463661030356312037$ ($|$diff$|=0$); $\pi_{1.5}=\pi_3=3.25976799305899509713550750411371829718207468$ ($3.6\times10^{-71}$); $\pi_{2.5}=3.18508354891312726234\ldots$ ($0$); $\pi_4=3.39693482362845842121\ldots$; $\pi_6=3.57268318383772505080\ldots$; $\pi_{10}=3.73396080922282829958\ldots$ Raw integral vs cosh: $2.2\times10^{-37}$ ($p=2$), $1.1\times10^{-50}$ ($p=1.5$). Cosh vs split form: $4.2\times10^{-81}$.

**Coefficients** (from exact $J$-integral arithmetic at 110 dps; first 11 cross-checked vs independent 150-digit direct fit, agreement $3.5\times10^{-20}$ for $N_0$):

| $k$ | $N_k$ (40 digits) | | $k$ | $N_k$ |
|---|---|---|---|---|
|0|−0.6931471805599453094172321214581765680755||5|−0.09559446392183336242739502251067627455892|
|1|0.2402265069591007123335512631633324858653||6|−0.9326682134840535164505582815863254270674|
|2|0.3557294080472350291649615278928845399454||7|−1.472571756898086547743541305319946712566|
|3|0.4363205189488704682613724388649403241718||8|0.1221389382135391530907648898935863780245|
|4|0.3496802380198935429771012951540838252324||9|4.700293562406214840901019828506766032978|

then $N_{10}=4.4249$, $N_{11}=-18.181$, $N_{12}=-30.262$, $N_{13}=149.34$, $N_{14}=299.17$, $N_{15}=-2073.5$, $N_{20}=-2.5866\times10^{6}$, $N_{25}=7.3762\times10^{10}$, $N_{30}=8.7559\times10^{14}$, $N_{35}=-2.4770\times10^{20}$, $N_{40}=-1.2405\times10^{25}$ (full 40-digit table in step_05 output).

**Growth:** $|N_k|^{1/k}=1.16,1.50,1.89,2.30,2.72,3.15,3.58,4.02,4.24$ at $k=10,14,18,22,26,30,34,38,40$ vs $k/(e\pi)=1.17,1.64,2.11,2.58,3.04,3.51,3.98,4.45,4.68$; odd-$k$ normalized $N_k\pi^{k-1}/(2(k-2)!)=1.72897,1.40548,1.28331,1.21873,1.17879,1.15161,1.13190,1.12393$ at $k=13,17,21,25,29,33,37,39$ → target $(-1)^{(k-1)/2}$ magnitude 1.

**Closed-form residuals:** $N_2$ (PSLQ, $1.1\times10^{-97}$), $N_3$ ($6.5\times10^{-97}$), $N_4$ ($2.9\times10^{-97}$; direct $1.9\times10^{-102}$), $N_5$ ($5.0\times10^{-102}$); $J_{3,0},J_{2,1},J_{4,0},J_{3,1},J_{2,2}$ PSLQ residuals $3.3\times10^{-101}$–$7.1\times10^{-102}$.

**Optimal truncation (asymptotic character):** $p=10$ ($h=0.1$): truncation errors $3.3\times10^{-6}$ ($k=3$), $2.0\times10^{-11}$ ($k=10$), $7.9\times10^{-17}$ ($k=31$, best) — least term at $k\approx\pi/h$, error $\sim e^{-\pi/h}$. At $p=1.1$: best error $1.1\times10^{-16}$ at $n=34$ ($\pi/h=34.6$).

---

## 4. Scripts (in working directory `D12_pi_p/`)

- **step_01.py** — Task 1: raw integral vs cosh form; $\pi_2=\pi$, $\pi_1=\pi_\infty=4$; duality table (60–70 digits); minimum at $p=2$.
- **step_02.py** — exact split $\pi_p=4\cdot2^{-1/p}+R$ (80-digit check); hypergeometric evaluation $g(h)=4\cdot2^{-h}$; $J(1,0),J(1,1),J(2,0)$ evaluations.
- **step_03.py** — closed forms $N_0..N_3$; confirmation against $\pi_p$ at $p=40,80,160$.
- **step_04.py** — exact series arithmetic: $[h^n]\Psi$ as rational polynomials in $(s,L,\log2)$ (uses $\log(1+e^{-x})=\log2-\frac x2+\log\cosh\frac x2$), saved to `/tmp/psi.pkl`.
- **step_05.py** — assembly of $N_0..N_{40}$ from 969 $J$-integrals at 110 dps; cross-check vs 150-digit direct fit.
- **step_06.py** — PSLQ identifications of $N_2..N_5$ and $J_{3,0},J_{2,1},J_{4,0},J_{3,1},J_{2,2}$ (95–100 dps); $\sqrt2$-absence test.
- **step_07.py** — growth/radius analysis; integrand-singularity asymptotics check by Cauchy integration; optimal truncation at $p=10$.
- **step_08.py** — $p=1$ expansion in $\eta=\frac1q-1$ and in $\eta=q-1$ ($M_n$ table, growth, verification); Euclidean-constant comparison.

(Note: files were executed in this session's environment where mpmath/sympy are installed; a bare `python3` without those packages will need them installed to re-run.)

## 5. Summary of labels

- **PROVED:** exact split $\pi_p=4\cdot2^{-1/p}+R(p)$; cosh representation; $N_0,N_1,N_2,N_3$ closed forms; exact formula for every $N_k$ (rational combination of alternating MZVs); $J_{1,k}=k!\,\eta(k+2)$, $J_{2,0}=\zeta(3)/4$.
- **PROVED modulo standard Euler sums** (all verified $\ge100$ digits): $N_4,N_5$ closed forms.
- **NUMERICAL:** duality/formula checks (40–80 digits); radius $=0$ with $|N_k|\asymp2(k-2)!/\pi^{k-1}$ and period-4 signs; $e^{-\pi/h}$ optimal-truncation scale; the $a_3..a_6$ pattern; $\sqrt2$-absence.
- **CONJECTURE:** nothing left unlabeled — the only unproven analytic claim is the rigorous large-order asymptotics (mechanism identified and numerically validated).
---

## Planner note, 2026-11-22 (after the referee pass REF_PIP_referee)

Three corrections to this file. Details and code: `../REF_PIP_referee/REFEREE.md`.

1. **All-orders remainder.** Section 2 proves the expansion of E(p) only through h^4 (h = 1/p). For N_4 and N_5 the same bound is needed to order h^6. It closes as follows. Write w = W(h, sh) with W(x, y) = e^{-xy/(1-x)} (1+e^{-y/(1-x)})^{-1-x}. W is analytic for |x| < 1/2, |y| < pi/2, so Cauchy estimates bound its Taylor remainder of total degree > m by C (h(1+s))^{m+1} for s <= h^{-1/2}; this is O(h^{m+1}) after integration against L = log(1+e^{-s}). The tails beyond h^{-1/2} and the expansion of e^{hL} - 1 are handled as at order 4.
2. **J_{2,2} is not used.** N_4 needs J_{2,1} and J_{3,0}; N_5 needs J_{3,1} and J_{4,0} (plus the lower ones). The five nontrivial sums are J_{2,0}, J_{2,1}, J_{3,0}, J_{3,1}, J_{4,0}.
3. **Growth law.** |N_k| ~ 2(k-2)!/pi^{k-1} holds only for odd k (ratio 1.28, 1.22, 1.18 at k = 21, 25, 29). For even k the ratio is about 11/k (0.51, 0.43, 0.375 at k = 22, 26, 30). Signs have period four. Divergence of the series is numerical, not proved.

Closed forms confirmed exactly (double shuffle + duality, independent code) and to 165 digits:
- N_4 = -log^5 2/120 + (pi^2/24)(1 - log 2 + log^2 2/2) + (zeta(3)/4)(1 - log 2) + pi^4/2880
- N_5 = log^6 2/720 + (pi^2/24)(1 - log 2 + log^2 2/2 - log^3 2/6) + (zeta(3)/4)(1 - log 2 + log^2 2/2) + (pi^4/2880)(1 - log 2) - 7 pi^4/1920 - zeta(3) pi^2/48 + zeta(5)/8
