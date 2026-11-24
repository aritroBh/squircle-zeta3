# SRF4 — Adversarial referee report

**Object.** `proof.md` in this folder: the claim that, for the area $S$ of $|x|^p+|y|^p+|z|^p=1$ expanded in
$\tau=1-1/p$, $a_k=[\tau^k]S/S(1)$, $z_0=(-\log2+i\pi)/2$, $R=|z_0|$:
$\limsup_k|a_k/k!|^{1/k}=1/R$; the Borel transform is analytic in $|z|<R$ and singular at $z_0,\bar z_0$; and
$a_k/k!=-\sqrt{3\pi/8}\,k^{-3/2}\operatorname{Re}z_0^{1/2-k}+o(k^{-3/2}R^{-k})$.

**Method.** I read every step of `proof.md` §§1–7 and redid each computation by hand. All numbers below come from
my own scripts in `referee/` (python3, mpmath 1.4.1, numpy). I did not run the author's scripts. I read the author's
coefficient file only after my own coefficients existed, and only to compare.

| script | what it does | output |
|---|---|---|
| `r1_representation.py` | $S$ from the definition (radial chart) vs $8(1-\tau)^2T(\tau)$ | `r1_out.txt` |
| `r2_lemmaZ.py` | all zeros of $1+e^{w\zeta}+e^{\zeta}$ in $\lvert\zeta\rvert<7$, 2007 values of $w$ | `r2_out.txt` |
| `r3_rouche_branch_hankel.py` | Lemma R (Rouché), Lemma B (branch), Lemma H (exact transfer) | `r3_out.txt` |
| `r5_coeffs.py K dps` | $a_k$ by an exact Dirichlet-moment expansion (own route) | `r5_coeffs_*.json` |
| `r5_law.py` | series vs quadrature; law, residual, fit; $d$-law | `r5_law_K40_out.txt`, `r5_law_K100_out.txt` |
| `r5_compare_author.py` | comparison with the author's `coeffs_K150_dps330.json` (run last) | `r5_compare_out.txt` |

---

## Item 1. The starting representation — HOLDS

*Derivation.* On the face in the positive octant, $z=(1-x^p-y^p)^{1/p}$ and
$dA=|\nabla F|/|F_z|\,dx\,dy$ with $F=\sum|x_i|^p$, so $dA=(\sum x_i^{2(p-1)})^{1/2}z^{1-p}\,dx\,dy$. Put
$t_i=x_i^p$, $dx_i=\frac1p t_i^{1/p-1}dt_i$. With $\frac1p=1-\tau$ and $\frac{p-1}p=\tau$:
$dA=(1-\tau)^2(t_1t_2t_3)^{-\tau}(\sum t_i^{2\tau})^{1/2}dt_1dt_2$. Eight octants give $S=8(1-\tau)^2T(\tau)$,
and $T(0)=\sqrt3/2$ gives $S(1)=4\sqrt3$. This is the proof's starting point.

*Lemma 1 (ordered log-ratio coordinates).* I rechecked the Jacobian $t_1t_2t_3$ and the identity
$\sum t_i^{2\tau}=t_1^{2\tau}P_w(\tau\ell)$, $f=t_1^{-2\tau}h_w(\tau\ell)$. Both are exact.

*Numbers* (independent route: radial chart $x=r(u)u$, $dA=r^3|\nabla F|/p\,d\sigma$):

| $\tau$ | $p$ | $S$ radial chart | $8(1-\tau)^2T$ | rel. diff. |
|---|---|---|---|---|
| 0.5 | 2 | 12.56637061435917 ($=4\pi$) | 12.56637061435144 | $6\cdot10^{-13}$ (quadrature limit) |
| 0.2 | 1.25 | 8.604826530243511 | 8.604826530243511 | $3\cdot10^{-20}$ |
| 0.02 | 1.0204 | 7.069313081596251 | 7.069313081596251 | $2\cdot10^{-21}$ |

---

## Item 2. Lemma Z and the upper bound — HOLDS

*Line by line.* (i) For $0<b<\pi$ both $\sin b>0$ and $\sin(wb)\ge0$, so $\operatorname{Im}E_w>0$; $b=0$ gives $E_w>0$.
Correct. (ii) The reduction to $\pi\le|b|<\rho^*$ is correct: $a=0$ with $|b|=\rho^*$ would need
$b\equiv\pm2\pi/3\pmod{2\pi}$, and $\rho^*=3.2171$ is not. In the case $a\ge0$: $|1+e^\zeta|=e^{wa}\ge1$ gives
$e^a\ge2\cos\beta$, so $a\ge\log(2\cos\beta)\ge0$ and $|\zeta|^2\ge g(\beta)$. The bound
$g'\ge2\pi-2\log2\tan\beta>0$ is right (it uses $\log(2\cos\beta)\le\log2$, $\tan\beta\ge0$). Equality forces
$\beta=0$, $a=\log 2$, and then $\operatorname{Im}E_w=e^{wa}\sin(w\pi)=0$ forces $w\in\{0,1\}$. The case $a<0$ follows from
$E_{1-w}(-\zeta)=e^{-\zeta}E_w(\zeta)$. I found no gap.

*Numbers* (`r2_out.txt`; 2007 values of $w$, including $10^{-6}$ and $1-10^{-6}$; Newton count = argument-principle
count at every $w$):
- $\min|\zeta|/\rho^*=1$ only at $w=0$ and $w=1$. On $0<w<1$ the minimum is $1+5.0\cdot10^{-7}$, at $w=10^{-6}$.
- Near the ends, $\min|\zeta|/\rho^*\approx1+w/2$ ($1.0005003$ at $w=10^{-3}$), which matches $z^{*\prime}(0)=z_0/2$.
- Profile symmetric under $w\mapsto1-w$ (to $2\cdot10^{-11}$); $1.302$ at $w=\frac12$.
- $\min|\operatorname{Im}\zeta|/\pi=1$; the inequality $a\ge\log(2\cos\beta)$ holds at all 78 zeros in the band.

*Does it give "no singularity nearer than $|z_0|$"?* Yes. The uniformity question disappears: Lemma 1 is exact,
so the zero-free disc of $f(\cdot;t)$ in $\tau$ is exactly $|\tau|<R/\ell(t)$ for every $t$, with no approximation in
$w$. Lemma U's Cauchy bound uses only $|t_1^{-2\tau}|\le3^{2|\tau|}$ and $|h_w(z)|\le\sqrt3e^{2R}$ for $|z|<R$, so the
constant does not depend on $w$ or $t$. The integration over $t$ is Lemma D (differentiation under the integral; I
checked the majorant $\sqrt3\,t_{\min}^{-1/2}(16+2\log\frac1{t_{\min}})^k$ and $\Pr(t_{\min}<\eta)\le6\eta$) plus
$\mathbb E[\ell^k]\le\mathbb E[\log^k(1/t_{\min})]\le6\,k!$. Together they give $|T_k|/k!\le C(R-\varepsilon)^{-k}$, so
the Borel transform has radius $\ge R$. The corners need no separate treatment: in region II the measure factor
$e^{-\lambda}\le e^{-w_0\ell}$ gives radius $\ge(1+w_0)(R-0.1)=1.795>R$. I checked this bound.

That no point of the circle other than $z_0,\bar z_0$ is singular is Proposition S(b). I checked it: on the compact
cut set, $\operatorname{Re}(z_1/\zeta)=1$ forces $\zeta=z_1$ with $|\zeta|=R$, hence $v=0$, $w=0$ and $\zeta=z_0$ (by
Lemma Z's equality case). That is exactly where "equality only at $w\in\{0,1\}$" is needed, and $w=1$ lies in
region II, which is analytic in a larger disc. Correct.

---

## Item 3. Lemma R, Lemma B, Lemma H — HOLDS

*Lemma R.* The zeros of $P_0=2+e^{-2z}$ and the moduli $R$ and $4.725$ are right. I recomputed the hand bound:
window $\varphi\in(1.39663,1.88935)$, $\sin\varphi>0.949690$, $|\psi|>3.798761\in[\pi,2\pi]$,
$\cos\psi>-0.791725$, $|P_0|\ge\sqrt{4-(2\cdot0.791725)^2}=1.221755$. The other two cases give $\ge1.5$. Rouché:
$e^{4\cdot0.19}-1=1.138276<1.221755$ (margin $0.083$). The true $\min_{|z|=2}|P_0|=1.495078$ (Lipschitz-certified
$\ge1.494735$). The largest ratio $|P_w-P_0|/|P_0|$ on the circle for $w\le0.19$ is $0.519$. The argument-principle
count is $2$ at all 381 values of $w\in[0,0.19]$. $|z^*|/R$ increases, to $1.108927$ at $0.19$. The finite
difference gives $z^{*\prime}(0)=-0.17324+0.78547i$ against $z_0/2=-0.17329+0.78540i$. And $v_1\in[0.1212,0.2433]$,
so $1+v_1\le1.25$ holds.

*Lemma B (sign of the root).* $D_0(sz_0)=(2+2^se^{-i\pi s})/(1-s)$ has $\operatorname{Im}<0$ on $(0,1]$, so the continued
root is the principal one. $(-z_0)^{1/2}=-iz_0^{1/2}$ and $e^{z_0}=i/\sqrt2$ give $\Psi^h_0(z_0)=\sqrt2z_0^{1/2}$. Correct.
Numerically, the root tracked along $[0,z^*]$ equals the principal root at $w=0,0.05,0.1,0.19$, with the argument
staying in $[-1.69,0]$. $\Psi^h_0(z_0)=1.123388474917+1.398266371667i=\sqrt2z_0^{1/2}$ (12 digits).

*Lemma H (exact remainder).* I rederived the keyhole decomposition. The positively oriented boundary of the
slit disc traverses the upper bank outward, where $(1-\zeta/\zeta^*)^{1/2}=-i\sqrt v$. This gives
$-\frac1\pi\zeta^{*-k}\int\sqrt v\,\Psi(\zeta^*(1+v))(1+v)^{-k-1}dv$ per cut. The test case agrees:
$-\frac1\pi\mathrm B(\frac32,k-\frac12)=[x^k](1-x)^{1/2}$ (checked at $k=2,5,9$). The bounds on $e_1$ (secant bound,
$\mathrm B(\frac52,k-\frac32)$) and $e_2$ hold. For $e_2$ the stated $(1+v_1)^{1-k}\frac\pi2$ is generous: the true
bound is $(1+v_1)^{1/2-k}/(k-\frac12)$. The circle term $M_F\rho_1^{-k}$ holds. The constants $M_h$, $L_h$, $M_\Psi$, $L$
are finite by joint continuity on the compact set $K_{\rm cut}$. The proof needs only that they are finite, not
their values. Numerically, $(\text{coef}-\text{lead})\cdot|z^*|^kk^{5/2}$ stays $O(1)$ for $k=20..80$ at
$(w,\ell)=(0,40),(0.05,3),(0.19,1)$. At $(0.19,1)$ the $\rho_1^{-k}$ term dominates for small $k$, as the bound allows.

*Lemma L and Theorem G.* I rechecked the bookkeeping. $T^{\rm III}_k$ gives
$M_k=-\frac{12}\pi\mathrm B\operatorname{Re}[z_0^{-k}I_k]$. Dominated convergence applies with majorant $M_\Psi$, since
$|z^*|\ge R$ gives $|(z^*/z_0)^{-k}|\le1$. Also $kg(\lambda/\ell)\to\lambda/2$, because $g'(0)=z^{*\prime}(0)/z_0=\frac12$ is
real, so the imaginary part is $o(1)$. The profile integral is $\int_{1/2}^1\sqrt{s(1-s)}ds=\pi/16$. The constants
check: $\frac{12}\pi\cdot\frac{\sqrt\pi}2\cdot\frac{\sqrt2\pi}{16}=\frac{3\sqrt{2\pi}}8$, and
$\frac2{\sqrt3}\cdot\frac{3\sqrt{2\pi}}8=\sqrt{3\pi/8}$.

---

## Item 4. The sharper law and the limsup — HOLDS WITH FIX (cosmetic)

Theorem G is proved with an **additive** error $o(k^{-3/2}R^{-k})$. Every error term is bounded in absolute value
on the scale $k^{-3/2}R^{-k}$ or smaller (regions I, II: $O(k!c^{-k})$ with $c>R$; Lemma H remainder:
$O(k^{-5/2}R^{-k})$; Lemma L: $o(1)$ inside $\operatorname{Re}[z_0^{-k}(\cdot)]$, which is $o(R^{-k})$). Nothing is
divided by $\operatorname{Re}z_0^{1/2-k}$. So the statement holds even where the cosine is small.

The limsup does follow. Corollary C's two-point lemma is correct: $\alpha+\beta\ge\operatorname{dist}(\theta,\pi\mathbb Z)=\pi-\theta$,
so $\max(|\cos x|,|\cos(x+\theta)|)\ge\sin\frac{\pi-\theta}2=0.62632$. So along at least one of every two consecutive
$k$, $|a_k/k!|\ge(c+o(1))k^{-3/2}R^{-k}$ with $c>0$, and $|a_k/k!|^{1/k}\to1/R$ along that subsequence.

**Fix (wording only).** Corollary C says "among any two consecutive $k$ one has
$|a_k/k!|\ge(\ldots)k^{-3/2}R^{-k}$" with the *same* $k$ on the right for both members. For the member $k+1$ the
bound is $(\ldots)(k+1)^{-3/2}R^{-k-1}$, which is smaller by a factor $R$. Replace it with: "for each $k$, one of
$j\in\{k,k+1\}$ has $|a_j/j!|\ge(0.6263\sqrt{3\pi/8}R^{1/2}+o(1))\,j^{-3/2}R^{-j}$." The conclusion is unchanged.

---

## Item 5. Independent numerics — HOLDS

*Route.* I did not use the author's Gamma-factorization engine. I expanded
$f=e^{-\tau\sum L_i}\sum_m\binom{1/2}m3^{1/2-m}X^m$ with $X^m=(y_1+y_2+y_3-3)^m$ and $y_i=t_i^{2\tau}$. Each monomial
integrates exactly by the Dirichlet formula
$\int_\Delta\prod t_i^{a_i}=\prod\Gamma(1+a_i)/\Gamma(3+\sum a_i)$, with $a_i=\tau(2j_i-1)$, and the $\Gamma$ series come
from $\log\Gamma(1+x)=-\gamma x+\sum_{n\ge2}(-1)^n\zeta(n)x^n/n$. The result is exact up to working precision.

*Validation.*
- $T_0=\sqrt3/2$ and $a_1=1$.
- $K=40$ at 120 and 170 digits agree.
- The truncated series matches direct quadrature of $S/S(1)$: at $\tau=0.02$ the difference is $2.6\cdot10^{-26}$,
  at the quadrature's precision; at $\tau=0.03$ it is $1.1\cdot10^{-24}$, the size of the next term.
- After that, comparison with the author's `coeffs_K150_dps330.json` (different engine): agreement to all
  110 compared digits for $k\le40$.

*Law, $k\le40$* (`r5_law_K40_out.txt`; residual $=(a_k/k!-\text{law})/(k^{-3/2}R^{-k})$):

| $k$ | $\lvert a_k/k!\rvert^{1/k}$ | $a_k/\text{law}$ | residual | $k\cdot$residual | $\lvert\cos((k-\frac12)\theta)\rvert$ |
|---|---|---|---|---|---|
| 10 | 0.418018 | 1.5018 | 0.1997 | 2.00 | 0.289 |
| 20 | 0.503874 | 1.0207 | 0.0271 | 0.54 | 0.953 |
| 25 | 0.518443 | 0.9852 | 0.0200 | 0.50 | 0.984 |
| 30 | 0.525096 | 0.9560 | −0.0477 | −1.43 | 0.789 |
| 33 | 0.491671 | **5.4387** | −0.0672 | −2.22 | 0.011 |
| 40 | 0.515374 | 1.6476 | −0.0550 | −2.20 | 0.062 |

The residual decays like $1/k$ while the ratio jumps (5.44 at $k=33$, the author's figure). This confirms that the
correct form is additive, not "$\sim$". Fit $k=20..40$ with three $1/k$-orders: amplitude $1.376725$ (theory
$\sqrt{3\pi R/8}=1.376611$) and phase $2.247541$ (theory $\pi-\theta/2=2.247616$). The sign and phase of Lemma B
are confirmed. The $d$-law residual is $0.0082$ at $n=40$ and decays.

$|a_k/k!|^{1/k}$ is still well below $1/R=0.62167$ at $k=40$, as expected from the $k^{-3/2}$ factor. With 40
coefficients this only supports the limsup; it cannot show it.

*Law, $k\le100$* (`r5_law_K100_out.txt`, 300 digits; $T_0-\sqrt3/2=3\cdot10^{-276}$; agrees with the author's file to
all 110 compared digits):

| $k$ | $\lvert a_k/k!\rvert^{1/k}$ | $a_k/\text{law}$ | residual | $k\cdot$residual |
|---|---|---|---|---|
| 60 | 0.563163 | 0.9877 | 0.0154 | 0.93 |
| 80 | 0.572680 | 1.0192 | 0.0189 | 1.51 |
| 91 | 0.531231 | **0.0218** | 0.0238 | 2.17 |
| 100 | 0.576369 | 0.9630 | −0.0199 | −1.99 |

The ratio at $k=91$ is $0.022$, the author's figure, while the residual stays small. Fit $k=20..100$ with three
$1/k$-orders: amplitude $1.376612$ against $1.3766111$ (6 digits), phase $2.247597$ against $2.2476162$ (5 digits).
Two-consecutive floor, each member on its own scale: $\min_{k\ge10}\max(\cdot)=0.8638$ and $\min_{k\ge60}=0.9108$, both
above the asymptotic floor $0.8622$ of Corollary C. The $d$-law residual is $0.0098$ at $n=100$, with $n\cdot$residual
bounded ($\approx1$–$1.6$). The ratio $R|a_k/k!|^{1/k}$ is $0.927$ at $k=100$; the remaining gap is the
$k^{-3/2}$ factor, $\approx e^{-1.5\log k/k}$.

*Not checked by me:* the $|B'(sz_0)|\sqrt{1-s}$ values in Proposition S(c) (illustrative, not used in any proof);
the paper's "120 coefficients … 9–10 digits" (I fitted only up to $k=100$).

---

## Item 6. The paper — wording stronger than what is proved, or now out of date

All phrases are on `paper/ieee/main.tex` line 1925 (one paragraph), plus Table II lines 2124–2125. The 15 lines
after 1925 are another topic (cuboctahedron) and need no change.

**Over-claims.**

1. Line 1925: "$a_k/k!\sim-\sqrt{3\pi/8}\,k^{-3/2}\operatorname{Re}z_0^{1/2-k}$". As a ratio this is false (ratio $5.44$ at
   $k=33$). Replace with: "$a_k/k!=-\sqrt{3\pi/8}\,k^{-3/2}\operatorname{Re}z_0^{1/2-k}+o(k^{-3/2}|z_0|^{-k})$. The error is
   additive. The ratio to the main term does not tend to $1$, because $\operatorname{Re}z_0^{1/2-k}$ comes close to $0$ for some $k$."
2. Line 1925: "the Borel transform in $\tau$ has square-root branch points at $z_0$ and $\bar z_0$ and none closer to the
   origin." Only the weak form is proved. Replace with: "the Borel transform in $\tau$ is analytic in $|z|<|z_0|$. On the
   circle $|z|=|z_0|$ it is singular only at $z_0$ and $\bar z_0$. Near them its derivative grows like
   $|1-z/z_0|^{-1/2}$, as at a square-root branch point. We have not proved that they are branch points in the strict
   sense."
3. Line 1925: "A Hankel transfer, which reads the coefficients off the nearest branch points, then gives the law." The
   proof does not read coefficients off branch points of the Borel transform. It applies the transfer to each point of
   the face. Replace with: "We get the law from a Hankel transfer at each point of the face, with error bounds that do not
   depend on the point, and then integrate over the face."
4. Line 1925: "Numerically, 120 coefficients confirm both laws to 9--10 digits". The law itself agrees only to $O(1/k)$
   on the scale $k^{-3/2}|z_0|^{-k}$. Only the fitted amplitude and phase can agree to many digits. Replace with:
   "Numerically, with 120 coefficients, the fitted amplitude and phase agree with the law to 9--10 digits." (I did not
   check the 9–10 digits. SRF4's own fit gives 5–7.)
5. Table II line 2125, column 3: "120 coefficients: law to 9--10 digits". Replace with: "120 coefficients: fitted
   amplitude and phase to 9--10 digits".

**Now under-stated or out of date (the proof is stronger than the text).**

6. Line 1925: "The sharper law is numerical, because its constant and phase rest on a Hankel transfer that is only
   sketched:". Replace with: "We also prove a sharper law:".
7. Line 1925: "the coefficients in $d$ follow the same law with $-\sqrt{3\pi/4}\,\operatorname{Im}$ in place of
   $-\sqrt{3\pi/8}\,\operatorname{Re}$." Replace with: "The coefficients in $d$ satisfy the same law, with
   $-\sqrt{3\pi/4}\,\operatorname{Im}$ in place of $-\sqrt{3\pi/8}\,\operatorname{Re}$, also with an additive error."
8. Line 1925: "The proof follows Theorem~\ref{thm:zero}." The new proof does not follow that route. Replace with: "We
   use the log-spread $\ell=\log(t_{\max}/t_{\min})$ of the face coordinates. In it the zero-free disc is exactly
   $|\tau|<|z_0|/\ell$, so no splitting into regions is needed."
9. Line 1925: "As for Theorem~\ref{thm:zero}, two uniformity steps are only sketched; the full argument, with an
   independent referee's repairs, is in our notes." Replace with: "Every step is written in our notes, and an
   independent referee has checked it."
10. Line 1925: "and the power $k^{-3/2}$ is sharp." This is now proved. Replace with: "The power $k^{-3/2}$ is sharp. We
    prove that, of any two consecutive $k$, one has $|a_k/k!|\ge c\,k^{-3/2}|z_0|^{-k}$ with $c>0$."
11. Table II line 2125: "Octahedron-end growth law, constant and phase (Sec.~\ref{sec:surfaces}-C) & Numerical (proof
    sketched) &". Replace with: "Octahedron-end growth law, constant and phase, additive error
    (Sec.~\ref{sec:surfaces}-C) & Proved$^\dagger$ &".
12. Table II line 2124: "Octahedron-end radius $0$; Borel singularities $z_0,\bar z_0$ … & Proved$^\dagger$ & refereed".
    This row is correct as stated, since it does not say "branch points". Keep it, but "refereed" should now refer
    to SRF4 and this report, not to SRF3R.

Line 1925's "Then$^\dagger$ $\limsup_k|a_k/k!|^{1/k}=1/|z_0|=\ldots=0.621668\ldots$, so the series in $\tau$, and the one in
$d$, has radius zero" is proved as stated. "They come from the edges of the octahedron, through the equation
$1+2e^{2z}=0$" is correct ($P_0=0\iff1+2e^{2z}=0$, and $w\to0$ is the edge).

---

## Smaller remarks on `proof.md`

- §7, "What is not proved": stated honestly. Proposition S(c) shows that $z_0$ is singular (unbounded $B'$) with
  square-root-type coefficients. It does not show a local form $A+(1-z/z_0)^{1/2}C$.
- Lemma H, $e_2$: the bound is correct but loose (see Item 3). No change is needed.
- Status table: "PROVED (here)" for every listed item, except S′, is confirmed by this report.

## Verdict

**Yes. The claim $\limsup_k|a_k/k!|^{1/k}=1/|z_0|$ is proved** (upper bound: Lemma Z + Lemma U, with no gap; lower bound:
Lemmas R, B, H, L + Corollary C, with one cosmetic index fix), and so is the additive-error growth law. The only
statement still not proved is the strict "square-root branch point".
