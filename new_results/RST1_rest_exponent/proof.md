# RST1: the rest exponent $\sqrt{2\pi}$ in the decay law of $M_k$ — proof

For every integer **N = k+1 ≥ 100**,
\[
\boxed{\;|M_k-M_k^{\rm main}|\;\le\;1030\;N^{-1/2}\,e^{-\sqrt{2\pi N}}.\;}\tag{1}
\]
In particular $|M_k-M_k^{\rm main}|\le 1030\,e^{-\sqrt{2\pi N}}$, i.e. the rest term of
Observation `obs:Mk` is $O(e^{-\sqrt{2\pi N}})$ with explicit constant $1030$.
Since $\sqrt{2\pi}=2.5066282746\ldots>2$, this improves CX6's
$|M_k-M_k^{\rm main}|\le 1.06\times10^{7}e^{-2\sqrt N}$ to the **true rest exponent**,
the last numerical piece of the decay law. The certified raw constant is
$C_{\rm raw}=1025.8898389\ldots$; we rounded up to $C=1030$ (margin $4.11$, certified).

**Which singularity produces the rest (item 1).** In CX6's representation
$M_k=-J_k+\frac2\pi\operatorname{Im}\int_\Gamma H(s)K(s)P_k(s)\,ds$, the factor
$E(s)=e^{i\pi s/2}$ is the image of the branch point $z=i\pi/2$ of the Borel kernel
($s=1$ pole of $K$). The main term pairs $E$ with the pole of
$Q_k=s^k/(1+s)^{k+1}$ at $s=-1$ (saddle $r=\sqrt{2N/\pi}$, decay $e^{-\sqrt{\pi N}}$).
The **rest is produced by pairing the same image $E$ with the next kernel pole at
$s=-2$** — the $m=2$ weight term (DCYA1 §2) — through the explicit piece
\[
A_k^{[2]}(s)=-\frac{2s^{k-1}(s+1)}{(s+2)^{k+1}}
\]
of CX6's correction kernel $A_k=P_k-Q_k$: saddle $r=2\sqrt{N/\pi}$, decay
$e^{-\sqrt{2\pi N}}$. The next branch points ($s=\pm3$, i.e. $z=\pm3i\pi/2$, via
$K-E\sim-E^3$) and the $s=-3$ pole ($m=3$) both give only $e^{-\sqrt{3\pi N}}$
($\sqrt{3\pi}=3.0699>\sqrt{2\pi}=2.5066$), certified subdominant in §7.
The leading asymptotics, derived from the $q=2$ saddle in §7:
\[
M_k-M_k^{\rm main}\sim R_k:=-\frac4\pi\operatorname{Im}\int_{\rm ray}
H(s)\,e^{i\pi s/2}\,\frac{2s^{k-1}(s+1)}{(s+2)^{k+1}}\,ds
\;\sim\;\frac{\sqrt\pi}{N}\,e^{-\sqrt{2\pi N}}\cos\sqrt{2\pi N}.\tag{2}
\]

## Lemma ledger

| # | Lemma | Status | Where |
|---|---|---|---|
| L0 | CX6 inputs: exact decomposition (3) below; $|A_k|$ bound (4); $|H|\le2$, $\sqrt{\pi/2}\,x^{-3/2}$, $|K|\le2e^{-ar}$, $|K-E|\le2e^{-3ar}$ on $\Gamma$; CX6 (17),(18) | **PROVED** (CX6, refereed; reused as input) | §1 |
| R1 | Modulus-exponent bound: $\frac{\rho^j}{|s+q|^{j+1}}\le\frac{\sqrt2}{q}\bigl(\frac{\rho}{\rho+q/\sqrt2}\bigr)^j$ and $-ar+j\log\frac{\rho}{\rho+q/\sqrt2}\le\tilde E(r)=-ar-\frac{\sqrt2\,jq}{2r+1+q/\sqrt2}$ | **PROVED (by hand)** | §2 |
| R2 | Gaussian saddle lemma $F_{H,\lambda}(j,q)$ (exact form, condition $r_\lambda>c_0$) | **PROVED (by hand)** + **CERTIFIED** numeric vs arb quadrature (ratio $\le0.005$) | §3; `step_18b.py` |
| R3 | Universal box lemma $F^{\rm box}_{H,\lambda}(j,q)$ | **PROVED (by hand)** + **CERTIFIED** numeric | §3; `step_18b.py` |
| R4 | $T_1$ assembly: $\frac{2}{\pi}|T_1|\le\frac8\pi[\,F_H(k,2)+\sum_i\frac{2^i}{i!}F_H(k-i,i+1)\,]$ split into $q{=}2$ / mid $[2,12]$ / mid2 $[13,k/2]$ (convexity) / tail ($>k/2$) | **PROVED (by hand)** + **CERTIFIED** constants | §4 |
| R5 | Convexity: $Q_i=e^{\frac{a(i+1)}{2\sqrt2}}e^{-\sqrt{\pi(k-i)(i+1)}}$ convex in $i$, so $Q_i\le Q_2$ on $[2,k/2]$ and $Q_i\le Q_{13}$ on $[13,k/2]$ for $k\ge99$ | **PROVED (by hand)** + **CERTIFIED** (grid; $h(99)=12.78>0$, $h(k)>0$ for $k\ge99$) | §4; `step_18b.py` |
| R6 | $T_2$: $\frac2\pi|T_2|\le\frac4\pi F_{H,3a}(k,1)$, decay $e^{-\sqrt{3\pi k}}$ | **PROVED** | §5 |
| R7 | $T_3$, $J_k$ elementary (CX6 (18); $|J_k|\le2^{-(N-1)}$) | **PROVED** | §5 |
| R8 | Monotonicity: every normalized piece $\mathcal P(k)$ is non-increasing for $k\ge99$ (log-derivative signs) | **PROVED (by hand)** + **CERTIFIED** (grid, 0 violations) | §6; `step_14.py` |
| R9 | Constant: $C_{\rm raw}=$ (sum of arb-evaluated pieces at $k=99$) $=1025.8898389\ldots<1030$ | **CERTIFIED (arb)** | `step_15.py` |
| R10 | Sanity: $|M_k-M_k^{\rm main}|<1030N^{-1/2}e^{-\sqrt{2\pi N}}$ at every $k\in[99,400]$ ($N\in[100,401]$), max ratio $0.215$ | **CERTIFIED (arb quadrature)** | `step_17.py` |
| R11 | Item 1: rest $=R_k(1+o(1))$; saddle amplitude $\to\sqrt\pi/N$, phase $\to\theta_2$; $q\ge3$ pieces subdominant | **CERTIFIED** (numerics) | §7; `step_19.py` |

## 1. Setup and CX6 inputs (L0)

Recall from CX6 (proof.md §§1–3, PROVED, reused verbatim):
\[
M_k-M_k^{\rm main}=-J_k+\frac2\pi\operatorname{Im}\{T_1+T_2-T_3\},\qquad
T_1=\int_\Gamma HKA_k\,ds,\quad T_2=\int_\Gamma H(K-E)Q_k\,ds,\quad
T_3=\int_0^{1/2}HEQ_k\,ds,\tag{3}
\]
with $\Gamma:\ s=\frac12+r\omega$, $\omega=e^{i\pi/4}$, $r\ge0$; $a=\frac{\pi}{2\sqrt2}$;
$\rho=|s|$, $x=\operatorname{Re}s=\frac12+\frac{r}{\sqrt2}$;
$|J_k|\le2^{-(N-1)}$;
\[
|A_k(s)|\le\frac{2\rho^k}{|s+2|^{k+1}}+2\sum_{i=1}^k\frac{2^i}{i!}\,
\frac{\rho^{k-i}}{|s+i+1|^{k-i+1}}\qquad\text{[CX6 (12)]}\tag{4}
\]
and on $\Gamma$ [CX6 (4),(13)]:
$|H(s)|\le2$, $|H(s)|\le\sqrt{\pi/2}\,x^{-3/2}$, $|K(s)|\le2e^{-ar}$, $|K(s)-E(s)|\le2e^{-3ar}$.
Also CX6 (18): $\frac2\pi|T_3|\le\frac{2}{\pi k}3^{-k}$.

Define, for $j\ge0$, $q\ge1$, $\lambda>0$,
\[
F_{H,\lambda}(j,q):=\int_0^\infty e^{-\lambda r}\,|H(s)|\,\frac{\rho^j}{|s+q|^{j+1}}\,dr,
\qquad F_H:=F_{H,a}.
\]

## 2. The modulus-exponent lemma (R1) — PROVED

For $s=\frac12+r\omega=\frac12+\frac{r}{\sqrt2}+i\frac{r}{\sqrt2}$:
$r\le\rho\le r+\frac12$ (since $\rho^2=r^2+\frac{r}{\sqrt2}+\frac14$), and
$|s+q|\ge\rho+\frac{q}{\sqrt2}$: indeed
$|s+q|^2-\rho^2=(\frac12+\frac{r}{\sqrt2}+q)^2-(\frac12+\frac{r}{\sqrt2})^2
=q+\sqrt2\,qr+q^2$, while $(\rho+\frac{q}{\sqrt2})^2=\rho^2+\sqrt2q\rho+\frac{q^2}{2}$;
their difference is
$q+\frac{q^2}{2}+\sqrt2 q(r-\rho)\ge q+\frac{q^2}{2}-\frac{q}{\sqrt2}
=q\bigl(1-\frac1{\sqrt2}\bigr)+\frac{q^2}{2}>0$
(using $r-\rho\ge-\frac12$). Hence
$\frac{\rho^j}{|s+q|^{j+1}}\le\frac1{\rho+q/\sqrt2}\bigl(\frac{\rho}{\rho+q/\sqrt2}\bigr)^j
\le\frac{\sqrt2}{q}\bigl(\frac{\rho}{\rho+q/\sqrt2}\bigr)^j$.
With the Padé bound $\log(1+u)\ge\frac{2u}{2+u}$ ($u\ge0$; equality at $u=0$, and
$\frac1{1+u}\ge\frac4{(2+u)^2}\iff u^2\ge0$),
\[
j\log\frac{\rho}{\rho+q/\sqrt2}
=-j\log\Bigl(1+\frac{q}{\sqrt2\rho}\Bigr)
\le-\frac{2jq}{2\sqrt2\rho+q}
\le-\frac{2jq}{2\sqrt2(r+\frac12)+q}
=-\frac{\sqrt2\,jq}{2r+1+\frac{q}{\sqrt2}},
\]
using $\rho\le r+\frac12$. So with $c_0:=1+\frac{q}{\sqrt2}$,
\[
-ar+j\log\frac{\rho}{\rho+q/\sqrt2}\le\tilde E(r):=-ar-\frac{\sqrt2\,jq}{2r+c_0}.
\tag{5}
\]

## 3. The two $F_H$ lemmas (R2, R3)

Write $\tilde E_\lambda(r)=-\lambda r-\frac{\sqrt2\,jq}{2r+c_0}$ (the rate-$\lambda$ form).
Set $r_\lambda:=\sqrt{\frac{jq}{\sqrt2\lambda}}$, $c_0=1+\frac{q}{\sqrt2}$,
$\Sigma:=\sqrt{2\sqrt2\,\lambda jq}$. By AM–GM, for all $r\ge0$,
\[
\lambda r+\frac{\sqrt2\,jq}{2r+c_0}
=\frac{\lambda}{2}(2r+c_0)+\frac{\sqrt2\,jq}{2r+c_0}-\frac{\lambda c_0}{2}
\ge2\sqrt{\frac{\lambda\sqrt2\,jq}{2}}-\frac{\lambda c_0}{2}
=\Sigma-\frac{\lambda c_0}{2}.\tag{6}
\]
Also $2\lambda r_\lambda=\Sigma$, and $\tilde E_\lambda$ is concave:
$\tilde E_\lambda''(r)=-\frac{8\sqrt2\,jq}{(2r+c_0)^3}<0$; its maximum is at
$r_m=\frac12(\sqrt{\frac{2\sqrt2\,jq}{\lambda}}-c_0)=r_\lambda-\frac{c_0}{2}$, with
$\tilde E_\lambda(r_m)=-\Sigma+\frac{\lambda c_0}{2}$.

**R3 (box lemma, universal).** Split $[0,2r_\lambda]$ and $[2r_\lambda,\infty)$.
On the first, $\tilde E_\lambda\le-\Sigma+\frac{\lambda c_0}{2}$ by (6); on the second,
drop the negative log term: $\tilde E_\lambda\le-\lambda r$. With $|H|\le2$ and R1,
\[
\boxed{F_{H,\lambda}(j,q)\le\frac{2\sqrt2}{q}\Bigl(2r_\lambda e^{\lambda c_0/2}
+\frac1\lambda\Bigr)e^{-\Sigma}.}\tag{7}
\]

**R2 (Gaussian lemma; condition $r_\lambda>c_0$).** Split
$[0,r_b]$, $[r_b,r_t]$, $[r_t,\infty)$ with $r_b=\frac{r_\lambda-c_0}{2}>0$, $r_t=2r_\lambda$.

*Left $[0,r_b]$:* $\tilde E_\lambda$ is concave and increasing on $[0,r_m]$; at $r_b$,
$2r_b+c_0=r_\lambda$, so
$\tilde E_\lambda(r_b)=-\frac{\lambda(r_\lambda-c_0)}{2}-\frac{\sqrt2\,jq}{r_\lambda}
=-\frac{\lambda r_\lambda}{2}+\frac{\lambda c_0}{2}-\Sigma
=-\frac54\Sigma+\frac{\lambda c_0}{2}$
(using $\frac{\lambda r_\lambda}{2}=\frac{\Sigma}{4}$ and $\frac{\sqrt2\,jq}{r_\lambda}=\Sigma$), and
$\lambda_L:=\tilde E_\lambda'(r_b)=-\lambda+\frac{2\sqrt2\,jq}{r_\lambda^2}
=-\lambda+2\sqrt2\,\lambda\cdot\sqrt2=3\lambda>0$.
Since $\tilde E_\lambda(r)\le\tilde E_\lambda(r_b)-\lambda_L(r_b-r)$ for $r\le r_b$
(concavity, $\tilde E_\lambda'\ge\lambda_L$ there), and $x^{-3/2}\le2^{3/2}$,
\[
\int_0^{r_b}e^{\tilde E_\lambda}x^{-3/2}dr\le\frac{2^{3/2}}{3\lambda}\,
e^{-\frac54\Sigma+\frac{\lambda c_0}{2}}.\tag{8}
\]

*Middle $[r_b,r_t]$:* here $x=\frac12+\frac{r}{\sqrt2}\ge\frac{r_b}{\sqrt2}
=\frac{r_\lambda-c_0}{2\sqrt2}=:x_b$.
On $[0,r_t]\supset[r_b,r_t]$, $\tilde E_\lambda''\le-c_2$ with
$c_2=\frac{8\sqrt2\,jq}{(2r_t+c_0)^3}=\frac{8\sqrt2\,jq}{(4r_\lambda+c_0)^3}$.
By Taylor (for $r$ between $r_m$ and the interval, using $\tilde E_\lambda'(r_m)=0$ and
$\tilde E_\lambda''\le-c_2$):
$\tilde E_\lambda(r)\le\tilde E_\lambda(r_m)-\frac{c_2}{2}(r-r_m)^2$. Hence
\[
\int_{r_b}^{r_t}e^{\tilde E_\lambda}x^{-3/2}dr
\le x_b^{-3/2}e^{-\Sigma+\frac{\lambda c_0}{2}}\int_{r_b}^{r_t}e^{-\frac{c_2}{2}(r-r_m)^2}dr
\le x_b^{-3/2}\sqrt{\frac{2\pi}{c_2}}\,e^{-\Sigma+\frac{\lambda c_0}{2}}.\tag{9}
\]

*Right $[r_t,\infty)$:* drop the log term: with $x\ge x_t:=\frac12+\frac{2r_\lambda}{\sqrt2}\ge\sqrt2\,r_\lambda$,
\[
\int_{r_t}^{\infty}e^{\tilde E_\lambda}x^{-3/2}dr
\le x_t^{-3/2}\int_{r_t}^\infty e^{-\lambda r}dr
\le\frac{(\sqrt2\,r_\lambda)^{-3/2}}{\lambda}e^{-2\lambda r_\lambda}
=\frac{(\sqrt2\,r_\lambda)^{-3/2}}{\lambda}e^{-\Sigma}.\tag{10}
\]

Assembling with $|H|\le\sqrt{\pi/2}\,x^{-3/2}$ and R1
($F_{H,\lambda}\le\frac{\sqrt2}{q}\int e^{\tilde E_\lambda}x^{-3/2}dr\cdot\sqrt{\pi/2}$;
note $\sqrt{\pi/2}\cdot\sqrt2=\sqrt\pi$):
\[
\boxed{F_{H,\lambda}(j,q)\le\frac{\sqrt\pi}{q}
\Bigl[\frac{2^{3/2}}{3\lambda}e^{-\Sigma/4}
+\Bigl(\frac{r_\lambda-c_0}{2\sqrt2}\Bigr)^{-3/2}\sqrt{\frac{2\pi(4r_\lambda+c_0)^3}{8\sqrt2\,jq}}
+\frac{(\sqrt2\,r_\lambda)^{-3/2}}{\lambda}\Bigr]
e^{\frac{\lambda c_0}{2}}e^{-\Sigma}}\tag{11}
\]
valid for $r_\lambda>c_0$. *Remark:* in (9) we used
$\sqrt{2\pi/c_2}=\sqrt{2\pi(4r_\lambda+c_0)^3/(8\sqrt2\,jq)}$; in the theorem we keep
this exact form (no relaxation), so the mid piece is
$\propto\frac1{r_\lambda}\bigl(\frac{4r_\lambda+c_0}{r_\lambda-c_0}\bigr)^{3/2}$.

**Numerical certification (R2/R3).** `step_18b.py` encloses the exact
$F_{H,\lambda}(j,q)$ by arb quadrature (saddle-aware segments $+$ proved tail) and
compares with (7),(11) for all $(j,q,\lambda)$ used in the proof and extras:
worst ratio exact/bound $=0.0050$ (box) and $0.0017$ (Gaussian). $\square$

## 4. The $T_1$ assembly (R4, R5)

From (3), (4), $|K|\le2e^{-ar}$ and the definitions,
\[
\frac2\pi|T_1|\le\frac8\pi\Bigl[F_H(k,2)+2F_H(k-1,2)
+\sum_{i=2}^{k}\frac{2^i}{i!}F_H(k-i,i+1)\Bigr].\tag{12}
\]
The first two terms are the $q=2$ piece ($i=0$ and $i=1$ both have $|s+2|$): they carry
the target decay $e^{-\sqrt{2\pi k}}$. Apply R2 ($r_\lambda>c_0$ holds: checked below).
The remaining sum is split:

*Mid $i\in[2,12]$:* R2 termwise (condition: worst case $i=12$, $j=k-12$,
$r_a=\sqrt{2(k-12)13/\pi}\ge26.8>c_0=1+13/\sqrt2=10.19$ at $k=99$). Each normalized term
is non-increasing (§6).

*Mid2 $i\in[13,k/2]$:* R3 (box) + convexity (R5). The exponent–prefactor product
$Q_i=e^{\frac{a(i+1)}{2\sqrt2}}e^{-\sqrt{\pi(k-i)(i+1)}}$ is **convex** in real $i$
($\ell''(i)=\sqrt\pi\bigl[\frac{1}{\sqrt u}+\frac{(k-2i-1)^2}{4u^{3/2}}\bigr]>0$ with
$u=(k-i)(i+1)$, after $\ell'(i)=\frac{a}{2\sqrt2}-\sqrt\pi\frac{k-2i-1}{2\sqrt u}$),
so $\max_{[13,k/2]}Q_i=\max(Q_{13},Q_{k/2})$. That $Q_{13}\ge Q_{k/2}$ for $k\ge99$
is $\ell(13)\ge\ell(k/2)$ with $\ell(i)=\frac{a(i+1)}{2\sqrt2}-\sqrt{\pi(k-i)(i+1)}$,
i.e. $h(k)\ge0$ with
$h(k)=\sqrt{\pi(k/2)(k/2+1)}-\sqrt{14\pi(k-13)}-\frac{a(k/2-13)}{2\sqrt2}$:
$h(99)=12.78>0$ (arb); $h(k)>0$ for $k\in[99,400]$ (direct arb check of
$Q_{13}\ge Q_{k/2}$) and for $k>400$ by
$h(k)\ge\frac{\sqrt\pi}{2}k-\sqrt{14\pi k}-\frac{ak}{4\sqrt2}
\ge0.69\,k-\sqrt{14\pi k}>0$ for $k\ge300$. With $(k-i)(i+1)\ge14(k-13)$ for $i\in[13,k/2]$ (downward parabola,
minimum at the left endpoint for $k\ge26$), $r_a\le\sqrt{k(k+1)/\pi}$ and
$\frac{2\sqrt2}{i+1}\le\frac{2\sqrt2}{14}$:
(the box prefactor is $e^{ac_0(i)/2}=e^{a/2}e^{\frac{a(i+1)}{2\sqrt2}}$, so the
prefactor–exponent product is $e^{a/2}Q_i\le e^{a/2}Q_{13}$ by R5, the factor $e^{a/2}$
being $i$-independent; the displayed $Q_{13}$ below already includes it):
\[
\sum_{i=13}^{k/2}\frac{2^i}{i!}F^{\rm box}_{H,a}(k-i,i+1)
\le\frac{2\sqrt2}{14}\Bigl(2\sqrt{\tfrac{k(k+1)}{\pi}}\,e^{\frac{a(1+14/\sqrt2)}{2}}
e^{-\sqrt{14\pi(k-13)}}+\frac1a e^{-\sqrt{14\pi(k-13)}}\Bigr)\sum_{i\ge13}\frac{2^i}{i!}.\tag{13}
\]

*Tail $i>k/2$:* R3, drop the exponent, CX6's (16) trick
($\sum_{i>k/2}\frac{x^i}{i!}\le e^{-k/2}e^{ex}$ for $e^{k/2}\le e^i$), with
$x'=2e^{a/(2\sqrt2)}$:
\[
\sum_{i>k/2}\frac{2^i}{i!}F^{\rm box}_{H,a}(k-i,i+1)
\le\frac{2\sqrt2}{k/2+1}\Bigl(2\sqrt{\tfrac{k(k+1)}{\pi}}\,
e^{\frac{a}{2\sqrt2}}e^{ex'}+\frac1a e^{2e}\Bigr)e^{-k/2}.\tag{14}
\]

## 5. $T_2$, $T_3$, $J_k$ (R6, R7)

$|K-E|\le2e^{-3ar}$ gives $\frac2\pi|T_2|\le\frac4\pi F_{H,3a}(k,1)$, and R2 applies
($r_{3a}=\sqrt{k/(\sqrt2\cdot3a)}\ge4.58>c_0=1+\frac1{\sqrt2}$ at $k=99$): decay
$e^{-\sqrt{3\pi k}}$, subdominant. CX6 (18): $\frac2\pi|T_3|\le\frac{2}{\pi k}3^{-k}$.
$|J_k|\le2^{-(N-1)}$.

## 6. Monotonicity (R8)

Each piece $\mathcal P$ of the assembled bound, normalized to
$t(k):=(k+1)^{-1/2}e^{-\sqrt{2\pi(k+1)}}$, is **non-increasing for $k\ge99$**.
The normalized pieces have log-derivatives
$\frac{d}{dk}\log\mathcal P=\frac{d}{dk}\log X_0+\frac{1}{2N}+\frac{\sqrt{2\pi}}{2\sqrt N}$
($N=k+1$), and in every case the negative part of
$\frac{d}{dk}\log X_0$ dominates $\frac1{2N}+\frac{\sqrt{2\pi}}{2\sqrt N}$:

* $q=2$ **mid** pieces ($i=0,1$): $X_0=\frac1{r_a}\bigl(\frac{4r_a+c_0}{r_a-c_0}\bigr)^{3/2}e^{-\sqrt{2\pi j}}$,
  $j\in\{k,k-1\}$: each factor $\frac1{r_a},\frac{4r_a+c_0}{r_a-c_0}$ is non-increasing in $k$,
  and the exponent part gives $\frac{\sqrt{2\pi}}{2\sqrt j}>\frac{\sqrt{2\pi}}{2\sqrt N}$,
  $\frac1{2j}>\frac1{2N}$ — so the log-derivative is negative termwise.
* $q=2$ **left** pieces: same, with an extra $-\frac{\sqrt{2\pi j}}{4}$ factor (more negative).
* $q=2$ **right** pieces: $X_0=r_a^{-3/2}e^{-\sqrt{2\pi j}}$: negative termwise
  ($\frac3{4k}>\frac1{2N}$ since $3N>2k$).
* **mid** ($i\in[2,12]$): exponent part $-\frac{\sqrt\pi(i+1)}{2\sqrt{(k-i)(i+1)}}$;
  $\frac{\sqrt\pi(i+1)}{2\sqrt{(k-i)(i+1)}}\ge\frac{\sqrt{3\pi}}{2\sqrt{k-2}}=0.1558\ldots$
  at $k=99$ vs $\frac1{2N}+\frac{\sqrt{2\pi}}{2\sqrt N}=0.1303\ldots$, and the ratio only
  improves for larger $k$; the R2-bracket factors are non-increasing in $r_a$ (hence in $k$).
* **mid2**: exponent part $-\frac{\sqrt{14\pi}}{2\sqrt{k-13}}=-0.3575\ldots$ at $k=99$
  dominates $\frac1{2N}+\frac{\sqrt{2\pi}}{2\sqrt N}+\frac{d}{dk}\log\bigl(\sqrt{\frac{k(k+1)}{\pi}}\bigr)$.
* **tail**: exponent $-\frac12$ dominates everything.
* **$T_2$**: exponent part $-\frac{\sqrt{3\pi}}{2\sqrt k}$; $\sqrt{3\pi}>\sqrt{2\pi}$.
* **$T_3,J_k$**: trivially super-exponential.

Grid certification (`step_14.py`): 0 violations of $\mathcal P(k+1)\le\mathcal P(k)$
on $k\in[99,1200)$ for each of $q2$, mid, $T2$, and on $[99,800)$ for mid2, tail, T3J.

## 7. Item 1: the producing singularity and the leading asymptotics (R11)

**The $q=2$ piece of $A_k$.** From CX6 (10), for large $|s|$ the $i=0$ (denominator) and
$i=1$ ($\beta$-tail) terms of $A_k$ are, in the saddle region $L\asymp\sqrt{\pi N}$
where $\beta(L)=\log(1+e^{-L})\approx e^{-L}$:
$A_k(s)\approx-\frac{s^k}{(s+2)^{k+1}}-\frac{s^{k-1}}{(s+2)^{k}}
=-\frac{2s^{k-1}(s+1)}{(s+2)^{k+1}}=A_k^{[2]}(s)$.
This is the pole of the weight kernel at $s=-2$ (the $m=2$ term of DCYA1's weight
expansion, an image of the *same* $\pm i\pi/2$ branch points — not the next branch
points $\pm3i\pi/2$). The resulting model for the rest is (2).

**Saddle analysis of $R_k$.** $R_k=-\frac4\pi\operatorname{Im}\int_{\rm ray}He^{\Phi(s)}ds$,
$\Phi(s)=\frac{i\pi s}{2}+(k-1)\log s+\log(s+1)-(k+1)\log(s+2)$. The saddle
$s_\star\approx2\sqrt{\frac{k-1}{\pi}}\omega$ solves
$\frac{i\pi}{2}s(s+1)(s+2)=s^2-2(k-1)s-2(k-1)$ (exact: $\Phi'(s_\star)=0$). Gaussian
evaluation gives, as $k\to\infty$:
\[
R_k\sim\operatorname{Im}\Bigl[-\frac4\pi\omega e^{\Phi(s_\star)}H(s_\star)
\sqrt{\frac{2\pi}{-\Phi''(s_\star)}}\Bigr]
\sim\frac{\sqrt\pi}{N}e^{-\sqrt{2\pi N}}\cos\sqrt{2\pi N},\tag{15}
\]
since $\operatorname{Re}\Phi(s_\star)=-\sqrt{2\pi(k-1)}-\frac12\log\frac{4(k-1)}{\pi}+O(N^{-1/2})$,
$H(s_\star)\sim\sqrt{\frac{\pi}{2}}s_\star^{-3/2}$,
$\sqrt{2\pi/(-\Phi'')}\sim\frac{2N^{1/4}}{\pi^{1/4}}$ (magnitude $\frac{\sqrt\pi}{N}e^{-\theta_2}$),
and the total phase $\arg(\cdot)-\theta_2\to\frac\pi2$ (so the Im gives $\cos\theta_2$).

**Certified numerics (`step_19.py`).**
- rest $=R_k\cdot(1+o(1))$: ratio $0.99971,1.00262,1.00004,0.99981,0.99999,0.99997,0.99999$
  at $k=100,150,200,250,300,350,400$ (exact $M_k$ from `mk_deg9.tsv`, quadrature $M_k^{\rm main}$).
- Saddle constants: $|g_k|Ne^{\theta_2}=1.9895,1.8651,1.8130,1.7904$ at
  $k=400,2000,10^4,5\times10^4$, $\to\sqrt\pi=1.77245$;
  $\arg(g_k)-\theta_2=1.51096,1.54405,1.55884,1.56545\to\pi/2$.
- Subdominance: $|R_k^{(q=3)}|/|R_k|=9.8\times10^{-4},4.8\times10^{-5},4.0\times10^{-6}$
  and $|R_k^{(K-E)}|/|R_k|=4.5\times10^{-3},2.5\times10^{-4},3.4\times10^{-5}$ at
  $k=100,200,400$ — both $O(e^{-(\sqrt{3\pi}-\sqrt{2\pi})\sqrt N})$.

## 8. The constant (R9) and the theorem

By R8, every piece of
$|M_k-M_k^{\rm main}|\le|J_k|+\frac2\pi(|T_1|+|T_2|+|T_3|)$
is non-increasing after normalization by $t(k)$ for $k\ge99$; hence
$|M_k-M_k^{\rm main}|\le C_{\rm raw}\,t(k)$ for all $k\ge99$ with
$C_{\rm raw}=$ the sum of the arb-evaluated pieces at $k=99$
(`step_15.py`, python-flint arb, 80 digits):

| piece | value at $k=99$ (arb, outward) |
|---|---|
| $q=2$ ($i=0,1$; (12) first two terms) | $1012.641802333060936174$ |
| mid ($i\in[2,12]$, R2 termwise) | $2.618057716968624948$ |
| mid2 ($i\in[13,k/2]$, (13)) | $5.649464712188422\times10^{-17}$ |
| tail ($i>k/2$, (14)) | $1.839745833717072\times10^{-5}$ |
| $T_2$ ($\frac4\pi F_{H,3a}(99,1)$) | $10.629960466490491331$ |
| $T_3+J_k$ | $1.213888250525348\times10^{-18}$ |
| **$C_{\rm raw}$** | $1025.889838913978389682$ |

$C_{\rm raw}<1030$ with margin $4.110$ (certified). This proves (1) for all $N\ge100$.

**Corollary (full $M_k$).** Combining (1) with DCY2's main-term theorem
($C_2=10^6$): the rest contributes
$1030N^{-1/2}e^{-\sqrt{2\pi N}}\le 667.357\,N^{-2}e^{-\sqrt{\pi N}}$ for $N\ge100$
(the function $1030N^{3/2}e^{-(\sqrt{2\pi}-\sqrt\pi)\sqrt N}$ is decreasing for
$\sqrt N>3/(\sqrt{2\pi}-\sqrt\pi)=4.086$, i.e. $N>16.7$; its value at $N=100$ is
$667.357$, arb). Hence with $\varpi_1$ of DCY2,
$$\Bigl|M_k+\frac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\frac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
\le C_3N^{-2}e^{-\theta},\qquad C_3 = 1\,000\,668\quad(N\ge100),$$
replacing DCY2's $C_3=C_2+1.85\times10^{10}$ (which used CX6's $e^{-2\sqrt N}$ rest).

## 9. Certified sanity table (R10)

`step_17.py` encloses $M_k^{\rm main}$ by arb quadrature (segments $[0,1,4,16,64,200]$
+ CX6's proved tail (20)) at 45 digits, subtracts the `mk_deg9.tsv` balls, and checks
(1) at **every** $k\in[99,400]$ ($N\in[100,401]$): **all pass**, max
$|M_k-M_k^{\rm main}|/t(k)=0.215$ at $k=106$. Printed also $|{\rm rest}|/{\rm lead}$
(${\rm lead}=\frac{\sqrt\pi}{N}e^{-\theta_2}$): values in $[0.09,1.18]$, oscillating as
$\cos(\theta_2+\varphi_N)$ — the predicted ratio. Sample:

| $k$ | $|M_k-M_k^{\rm main}|$ (arb) | $/$ $t(k)$ | $/$ lead |
|---|---|---|---|
| 99 | $1.64857408084\times10^{-13}$ | 0.127 | 0.716 |
| 150 | $4.403845828\times10^{-17}$ | 0.013 | 0.089 |
| 200 | $3.811407939\times10^{-18}$ | 0.147 | 1.174 |
| 250 | $1.884388565\times10^{-20}$ | 0.053 | 0.471 |
| 300 | $1.927871694\times10^{-22}$ | 0.026 | 0.252 |
| 350 | $1.220912698\times10^{-23}$ | 0.057 | 0.601 |
| 400 | $4.799125366\times10^{-25}$ | 0.061 | 0.684 |

## 10. Reproduction

Run from this directory with the project venv:
```sh
'<project>/.venv/bin/python' verify_rst1.py
```
It recomputes every certified constant (arb), validates R2/R3 against exact arb
quadrature, checks convexity and monotonicity on grids, rebuilds the constant $C$,
re-runs the sanity table $k\in[99,400]$, and re-checks the item-1 predictions.
Supplied inputs (`mk_deg9.tsv`, CX6's lemmas) remain read-only.
