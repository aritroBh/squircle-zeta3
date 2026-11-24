# D13s — Bug fixes B1/B2, proofs of the "harmless" items, and re-certified decay theorem

**Theorem (re-certified).** For all integers (equivalently reals) $N \ge 100$,
$$\Bigl|\,M_k^{\mathrm{main}} + \tfrac{\sqrt\pi}{N}\,e^{-\sqrt{\pi N}}\cos\sqrt{\pi N}\,\Bigr| \;\le\; C\,N^{-3/2}\,e^{-\sqrt{\pi N}},\qquad C = 254.2276 .$$

The constant is the max of three independently certified ranges:

| range | cover | certificate | sup $C$ |
|---|---|---|---|
| $[100,400]$ | 29 cells, ratio 1.05 | `fixed/step_18.py` | $194.1267251$ (worst cell $(373.346, 392.013)$) |
| $[400,10^{10}]$ | **85181 cells, ratio 1.0002** | `fixed/planner_sweep_ratio.py 1.0002` | $254.2275176$ (worst cell $[1890.7129, 1891.0910]$) |
| $[10^{10},\infty)$ | single analytic cell $y\in[0,y_J]$, $y_J=\sqrt{\pi/2\cdot10^{-10}}$ | `fixed/step_10_J192.py` | $47.0701250$ |

$$\boxed{C = \max(194.1267251,\;254.2275176,\;47.0701250) = 254.2275176 \;\;(\le 254.2276)}$$

Truth check (full contour, mpmath @50–60 dps): $|M^{\rm main}-{\rm lead}|/{\rm sc} = 0.602\ (N{=}100)$,
$1.065\ (N{=}400)$, $0.075\ (N{=}1891)$, $0.462\ (N{=}10^4)$, $0.473\ (N{=}10^5)$ — all $\ll C$.
Cover is contiguous ($400 = 400$, $10^{10} = 10^{10}$, $y_J \leftrightarrow N_J$) and **no cell
failed after the B2 fix** (29/29, 43/43, 85181/85181 accepted, 0 splits, final cell all asserts pass);
no subdivision was needed beyond the stated grids.

---

## Bug B1 (final cell `ecw`, off by $2^{-1/4}$) — FIXED

**Buggy line** (`step_10_J192.py:168`): `ecw = (rt2/pi).sqrt()*(-c_w/2).exp()` computes
$\sqrt{\sqrt2/\pi}\,e^{-c_w/2} = 0.4069446$; the tail identity
$e^{-c_2 d^2} = y\sqrt{2/\pi}\,e^{-c_w/2}(\dots)$ needs $\sqrt{2/\pi}\,e^{-c_w/2} = 0.4839414$
(the comment `# sqrt(2/pi)` shows the intent; `rt2` is $\sqrt2$, not $2$).
**Fix:** `ecw = (arb(2)/pi).sqrt()*(-c_w/2).exp()`. All of TM/TR/TLm grow by $2^{1/4}=1.189207$.

| quantity | before (buggy) | after (fixed) |
|---|---|---|
| TMreg | $0.10290245$ | $0.12237233$ |
| TM contribution to $C$ $= (2/\pi)(\pi/2)^{3/2}{\rm TMreg}$ | $0.1289691$ (**VIOLATED**: $<$ truth $0.1441235$ at $N=10^{10}$) | $0.1533710 > 0.1441235$ (margin $1.064$) |
| TR contribution | $0.29596$ (still $>0.11345$) | $0.35196 > 0.11331$ (margin $3.1$) |
| TLm contribution | $1.61775$ (still $>0.09084$) | $1.92384 > 0.09084$ (margin $21$) |
| final-cell $C$ | $49.5625979$ (**not a valid proof**) | $49.9490881$ with old breg; **$47.0701250$ with the H1 exact form** |

TM truth $=(2/\pi)\,{\rm TM}_{\rm erfc}/{\rm sc}$ is decreasing from the junction
($0.1441235$ at $10^{10}$, $0.1429748$ at $1.5{\cdot}10^{10}$, $0.1324781$ at $10^{12}$,
$0.1232537$ at $10^{14}$), so the junction $N=10^{10}$ is the worst case and the fixed bound
dominates everywhere on the final cell (`fixed/verify_TM_final.py`, `fixed/verify_TM_mono.py`).

## Bug B2 (all box cells: `A1m` does not enclose $A_1$) — FIXED

**Buggy line** (`step_38.py:165`): `A1m = iv(0, (z2/(1+z2)).upper())`.
The quantity required by the (exact) phase formula
$\varphi = \tfrac{\pi}{\sqrt2}B_h$, $B_h = -\tfrac{y\sigma}{2\rho} + \tfrac{2y^2\sigma-\sqrt2\rho y-2\rho A_1 y}{2y(2+\sqrt2\rho y)}$,
is $A_1 = \bigl(1-\tfrac{\arctan zz}{zz}\bigr)\Big/y$, and the coded box only encloses
$A_1' = 1-\arctan(zz)/zz$ (indeed $A_1' \le zz^2/(1+zz^2)$). For $y<1/3$ — every $N\ge100$ —
$A_1'/y > zz^2/(1+zz^2)$ is possible; confirmed: $\max\bigl(A_1'/y - zz^2/(1+zz^2)\bigr) = +0.0115$.

**Fix (proved enclosure).** $f(x) = 1 - \frac{\arctan x}{x}$ is nondecreasing on $[0,\infty)$:
$f'(x) = g(x)/x^2$ with $g(x) = \arctan x - \frac{x}{1+x^2}$, $g(0)=0$,
$g'(x) = \frac{2x^2}{(1+x^2)^2} \ge 0$ (sympy-printed). Since $f\ge0$ and $y>0$,
$$A_1 = \frac{f(zz)}{y} \in \Bigl[\frac{f(zz_{\rm lo})}{y_{\rm up}},\; \frac{f(zz_{\rm up})}{y_{\rm lo}}\Bigr],$$
with $zz = \tfrac{\sqrt2\rho y}{2+\sqrt2\rho y}$ boxed by arb over the cell (outer enclosure).
Implemented in `fixed/step_38.py`; containment of the true $A_1$ verified at 14112 sample
$(y,\sigma)$ points in 8 cells (worst violation $-9.9{\cdot}10^{-7}$, i.e. contained), and the full
$z$-box now contains $z_{\rm true}$ at 369/369 sampled $N\in[100,1.5{\cdot}10^9]$
(`fixed/verify_A1.py`, `fixed/verify_zbox.py`). Together with the exactness of
$\varphi=(\pi/\sqrt2)B_h[A_1'/y]$ (worst residual $1.9{\cdot}10^{-53}$, 400 pts @60 dps,
`fixed/verify_Bh.py`), this upgrades the review's numerical observation to a proof.

| range | buggy sup $C$ | fixed sup $C$ |
|---|---|---|
| $[100,400]$ (29 cells) | $194.3808879$ | $194.1267251$ |
| $[400,10^{10}]$ ratio 1.5 (43 cells) | $617117.5536$ | $617118.8072$ |
| $[400,10^{10}]$ ratio 1.0002 (85181 cells) | $254.5357$ | $254.2275176$ |

(The fine ratio-1.0002 cover is far tighter than the ratio-1.5 cover; it is the one quoted in the theorem.)

---

## H1 (D13h's `breg`) — exact regular form derived; final cell re-run with it

What the certificate needs: $|e^{i(\varphi+\pi/4)}-1| \le |\varphi+\pi/4| = \tfrac{\pi}{\sqrt2}y\,|B_{\rm true}| \le y\,d_6$
with $B_{\rm true} = \tfrac{\varphi+\pi/4}{(\pi/\sqrt2)y}$ on the saddle branch.

**Exact form (D13s).** On the saddle branch $\rho(t) = \sqrt{1+\tfrac{t^2}{1+\sqrt2 t}}$, $y=t/\rho$, $N=\pi/(2y^2)$
(saddle equation $\mathrm{Re}\,\Phi'(1/t)=0$ proved: $\rho^2(1+\sqrt2 t)=1+\sqrt2 t+t^2$), and using
$\arg(1+bx_0) = \arctan\tfrac{1}{1+\sqrt2 t}$ and $\tfrac{\pi}{4}-\arctan\tfrac{1}{1+\sqrt2 t}=\arctan(zz)$:
$$\varphi = \frac{a(1-2\rho)}{t} + N\arctan(zz),\qquad
B_{\rm true} = \frac{\rho\,B_1}{2+\sqrt2 t},\quad
B_1 = \tfrac12 + \sqrt2\,s_1 + s_1^2 - \rho^2 A_{1n},$$
$$s_1 = \frac{1-\rho}{t} = \frac{-t}{(1+\sqrt2 t)(1+\rho)},\qquad
A_{1n} = \frac{A_1'}{t^2} \in \Bigl[\frac{2}{3(2+\sqrt2 t)^2}-\frac{4t^2}{5(2+\sqrt2 t)^4},\; \frac{2}{3(2+\sqrt2 t)^2}\Bigr],$$
the $A_1'$ bracket from the alternating series $\arctan(z)/z = \sum(-1)^k z^{2k}/(2k+1)$ ($|z|\le1$).
Everything is manifestly regular at $t=0$ ($B_{\rm true}\to 1/6$).
**Verification:** worst relative residual vs direct $B_{\rm true}$ $= 8.5{\cdot}10^{-37}$ over 600 random
$t\in[10^{-12},1]$ @60 dps; the arb box over the final cell is $B_{\rm true}\in[0.16666\pm8.2{\cdot}10^{-6}]$,
containing the truth at 20000/20000 grid points (`fixed/verify_Btrue.py`).

**D13r's "exact form" is not exact.** D13r's step_31.py has an `IndentationError` (never ran), and its
step_32.py prints $B_{\rm true} - {\rm breg}_h + \tfrac1{12} \in [2.4{\cdot}10^{-9},\,0.044]$ — i.e.
${\rm breg}_h - \tfrac1{12}$ is *not* identically $B_{\rm true}$ (the result.md claim "residual 5.7e-44"
could not be reproduced). The review's conclusion (d6 valid) is nevertheless correct, and is now proved:
the coded box gives ${\rm breg\_abs} = 0.75 > B_{\rm abs} = 0.16666814$ (the printed
`assert breg_abs >= B_abs` passes), so the old $d_6$ was harmless. The fixed code replaces `breg` by the
exact form: $d_6 = \tfrac{\pi}{\sqrt2}B_{\rm abs} = 0.37024$ (was $1.66609$), ${\rm Dreg} = 8.56332$
(was $11.44228$), final-cell $C = 47.0701250$ ($49.9490881$ if the old breg is kept — also valid).

## H2 (TL0's constant) — proved bound substituted and re-run

The coded ${\rm TL0} = \sqrt2\,(2+\sqrt2)^{-N/2}e^{\theta}/(N-1)$ could not be derived by the review
(its own derivation gives $\sim2\times$). D13s proves and codes
$$I_0 := \int_0^1 e^{\mathrm{Re}\,\Phi(x)}|g(x)|\,dx \;\le\; (2+\sqrt2)^{-N/2}\frac{2}{N-2},
\qquad {\rm TL0} = (2+\sqrt2)^{-N/2}\,e^{\theta}\,\frac{2}{N-2},\quad \theta=\sqrt{\pi N}.$$
*Proof.* (i) $x^2|g(x)| = |s\,H_c(s)|\le1$ for $s = e^{i\pi/4}x$, $x\in(0,1]$: the beta-integral series
$H_c(s) = \tfrac12\int_0^1 u^{s/2-1}(1-u)^{1/2}du = \tfrac1s - \sum_{k\ge1}\frac{c_k}{s+2k}$
with $c_k = (-1)^{k+1}\binom{1/2}{k}>0$, $\sum c_k = 1$ gives
$sH_c(s) = 1-\sum c_k w_k$, $w_k = \tfrac{s}{s+2k}$. On the ray,
$2\,\mathrm{Re}\,w_k \ge |w_k|$ because
$[2(x+\sqrt2 k)]^2 - (x^2+2\sqrt2kx+4k^2) = 3x^2+6\sqrt2kx+4k^2 > 0$; since $|w_k|\le1$ also
$2\,\mathrm{Re}\,w_k \ge |w_k|^2$, hence
$|sH_c(s)|^2 = 1 - 2\sum c_k\mathrm{Re}\,w_k + \bigl|\sum c_k w_k\bigr|^2
\le 1 - \sum c_k|w_k|^2 + \sum c_k|w_k|^2 = 1$ (Jensen, $\sum c_k = 1$).
(ii) $1+\sqrt2x+x^2 \ge (2+\sqrt2)x$ (AM–GM $x+1/x\ge2$). Then
$I_0 \le (2+\sqrt2)^{-N/2}\int_0^1 x^{N/2-2}dx = (2+\sqrt2)^{-N/2}\frac{2}{N-2}$ ($N\ge100>4$). $\square$

Numerics: $\max x^2|g(x)| = 0.9999999995$ (grid), $2\mathrm{Re}(w)/|w| \ge 1.415$, truth/bound ratio
$\approx 0.27$ at $N=100,400,1000$ (`fixed/certify_H2.py`). Implemented in `fixed/step_38.py` and in the
final cell (${\rm TL0N}$ multiplied by $N^{3/2}$; decreasing in $N$ for $N\ge N_J$ since its log-derivative
$-\tfrac12\log(2+\sqrt2) + \tfrac{\sqrt\pi}{2\sqrt N} + \tfrac{1.5}{N} - \tfrac1{N-2}$ is negative at
$N_J=10^{10}$ and itself decreasing: its derivative is
$-\tfrac{\sqrt\pi}{4N^{3/2}} - \tfrac{1.5}{N^2} + \tfrac1{(N-2)^2} < 0$ for $N\ge 10^{10}$ as
$\sqrt\pi N^{1/2}(1-2/N)^2 > 4$). Contribution to $C$: $\sim10^{-18}$ (cells) / $10^{-2666376434}$ (final
cell) — $C$ unchanged, as predicted.

## H3 (mS) — restated proof (with a correction to D13r's margin formula)

With $t = \rho y$, $S(t) = \log(1+\sqrt2 t+t^2) - \sqrt2 t$, the exact identity (D13s, residual
$1.8{\cdot}10^{-31}$ over $t\in[10^{-10},10^{10}]$ @60 dps) is
$$\mathrm{Re}\,\Phi(x_0) + \theta = -\tfrac N2 S(t) \;-\; \frac{\pi\sqrt2}{4t}(\rho-1)^2
\qquad\text{(the correction term is NEGATIVE)}.$$
Write $-\tfrac N2S(t) = \tfrac N2\psi(t)$, $\psi(t) = \sqrt2 t - \log(1+\sqrt2t+t^2)$. Then
$\psi(0)=0$ and $\psi'(s) = \tfrac{\sqrt2 s^2}{1+\sqrt2 s+s^2}$, so for $t \le T$:
$\tfrac{\sqrt2}{3}t^3/(1+\sqrt2T+T^2) \le \psi(t) \le \tfrac{\sqrt2}{3}t^3$ — the h3-box
$[\,\tfrac{\sqrt2}{3}/(1{+}\sqrt2T{+}T^2),\; \tfrac{\sqrt2}{3}\,]$ brackets $\psi(t)/t^3$.
**Upper box (all $t>0$):** truth $= \tfrac N2\psi(t) - \tfrac{\pi\sqrt2}{4t}(\rho-1)^2 \le \tfrac N2\psi(t) \le \tfrac{N}{2}\tfrac{\sqrt2}{3}t^3$ (margin min $3.9{\cdot}10^{-21}$ over the scan).
**Lower box (certified domain):** with $H(t) := \psi(t) - \tfrac{\sqrt2}{3}\tfrac{t^3}{1+\sqrt2t+t^2} - \sqrt2\,t(\rho-1)^2/\rho^2$ ($=$ margin$/\tfrac N2$),
$H(0)=0$ and $H'(t) = t^3{\rm Hp3}(t)$ with the *manifestly regular*
$${\rm Hp3}(t) = \tfrac{\sqrt2}{3}\tfrac{\sqrt2+2t}{(1+\sqrt2t+t^2)^2} - \sqrt2 w^2 t - \tfrac{\sqrt2\,w\,t(2+\sqrt2t)}{q(1+\sqrt2t+t^2)^2},\quad
q = \sqrt{\tfrac{1+\sqrt2t}{1+\sqrt2t+t^2}},\; w = \tfrac{1}{(1+\sqrt2t+t^2)(1+q)}.$$
arb box evaluation certifies ${\rm Hp3}([0,0.13]) \ge 0.1538 > 0$ (`fixed/certify_H3.py`), hence $H\ge0$
on $[0,0.13]$; the certified domain $N\ge100$ has $t \le y(1+y^2)\big|_{y=\sqrt{\pi/200}} = 0.12731 < 0.13$.
($H$ stays positive up to $t^* = 4.5447$, i.e. down to $N = 0.288$; the lower box *does* fail for large $t$
— D13r's "for all $t>0$" is true only for the upper box, and D13r's margin formula
$\tfrac{\pi\sqrt2}{4t}[\tfrac{t^3(\sqrt2+t)}{3(1+\sqrt2t)} - (\rho-1)^2]$ is **not** a lower bound on the
true margin for small $t$: its bracket's first term $\sim\tfrac23 t^4$ exceeds the true
$G(t)\sim\tfrac16 t^4$. The containment conclusion on the certified domain is nevertheless true, as proved
above.) The arb-evaluated mS box then contains the truth (outer enclosure + monotonicity of the bounds),
and $z$-box containment is confirmed at 369/369 sampled $N$ with the fixed code.

## Smaller items (D13r list) also fixed in the re-run

* **TLfN missing $(\pi/2)^{3/2}$** (final cell): restored
  (`TLfN = G1h*(pi/2)**1.5*exp(-K_J/y_J)/(y_J^3*lam2)`); value $2.79{\cdot}10^{-19229} \to 5.50{\cdot}10^{-19229}$, no effect on $C$.
* D13r's step_30 print "Btrue inside breg-box: max violation $= 0.0833$" is, by its own criterion
  ("$<0$ or 0 $\Rightarrow$ contained"), a **non-containment**; harmless because only
  $|B_{\rm true}|\le{\rm breg\_abs}$ is needed — proved via H1 ($0.16667 < 0.75$).

---

## Lemma list and proof status

| lemma | content | status |
|---|---|---|
| L1 | $|\gamma'| \le \tfrac34(\tfrac1{t^2}+\tfrac1t)+\tfrac1x$, $|\gamma''| \le \tfrac38(\tfrac2{t^3}+\tfrac1{t^2})+\tfrac1{x^2}$ | **HOLDS** (D13r: proof + 1659 pts @60 dps) |
| L2 | $\lambda_R/\lambda_L/\lambda_2$ regular form $-\mathrm{Re}\,P_1 = a - \tfrac{\pi\rho}{2}G_h(\eta)$ | **HOLDS** (D13r: algebra + 900 pts, resid $9.5{\cdot}10^{-53}$) |
| L3 | $G_h$ decreasing in $\eta$; $|g|$ decreasing; $\beta$-bound $|g|\le H_r(x/\sqrt2)/x$ | **HOLDS** (D13r: proved, sympy) |
| L4 | $|\Phi_3|\le\Phi_{3b}$, $|\Phi_4|\le\Phi_{4b}$ + monotonicity | **HOLDS** (D13r: proved) |
| L5 | $\mathrm{Re}\,X_2 = \theta(M_1-L_2)$; $\log(1+u)/u$ decreasing | **HOLDS** (D13r) |
| L6 | Binet $|\mu(z)| \le 1/(12|z|\cos\arg z)$; real-$z$ $E$ bound; final-cell $\mathrm{Re}\,E\le0$ | **HOLDS** (D13r: proved) |
| L7 | final-cell factorisations, tail identity, $z$ = MV/$\tilde L$, leading-term identity | **HOLDS** (D13r), with D13s amendments: mS identity sign $(\rho-1)^2$ term is negative; h3-box containment proved on domain (H3); $A_1$ phase box replaced (B2) |
| L8 | assembly (coverage, W-decomposition, IBP, TL0) | **HOLDS** after D13s fixes B1 (ecw), H2 (TL0 constant), TLfN factor |
| D13s-A1 | $f(x) = 1-\arctan(x)/x$ nondecreasing on $[0,\infty)$; $A_1$ enclosure $[f(zz_{\rm lo})/y_{\rm up}, f(zz_{\rm up})/y_{\rm lo}]$ | **PROVED** (sympy $g'(x) = 2x^2/(1+x^2)^2 \ge 0$; 14112-pt containment) |
| D13s-Bh | $\varphi = \tfrac{\pi}{\sqrt2}B_h$ exact with $A_1 = A_1'/y$ | verified $1.9{\cdot}10^{-53}$ (400 pts @60 dps) |
| D13s-Btrue | exact regular form of $(\varphi+\pi/4)/((\pi/\sqrt2)y)$ | verified $8.5{\cdot}10^{-37}$ (600 pts); arb box containment 20000 pts |
| D13s-TL0 | $x^2|g(x)|\le1$ (series + $2\mathrm{Re}\,w_k\ge|w_k|$); $1+\sqrt2x+x^2\ge(2+\sqrt2)x$ | **PROVED** (above; truth ratio 0.27) |
| D13s-mS | identity + $\psi$ sandwich + $H\ge0$ via ${\rm Hp3}\ge0.1538$ on $[0,0.13]$ | **PROVED** (arb) |

## Evidence table (all numbers printed by code in this directory)

| check | file | result |
|---|---|---|
| B1: TM dominates truth at $N=10^{10}$ | `fixed/verify_TM_final.py` | bound $0.1533710$ vs truth $0.1441235$ (margin $1.064$) |
| B1: TM truth decreasing from junction | `fixed/verify_TM_mono.py` | $0.1441235 \searrow 0.1324781$ ($10^{10}\to10^{12}$) |
| B1: TR/TLm dominate | `fixed/verify_TM_final.py` | $0.35196>0.11331$ ($3.1\times$); $1.92384>0.09084$ ($21\times$) |
| B2: $A_1$ box containment | `fixed/verify_A1.py` | 14112 pts, worst violation $-9.9{\cdot}10^{-7}$ |
| B2: $z$-box containment (fixed code) | `fixed/verify_zbox.py` | 0 misses / 369 sampled $N\in[100,1.5{\cdot}10^9]$ |
| Bh exactness $A_1=A_1'/y$ | `fixed/verify_Bh.py` | worst residual $1.94{\cdot}10^{-53}$ @60 dps |
| H1: exact $B_{\rm true}$ form | `fixed/verify_Btrue.py` | worst rel resid $8.46{\cdot}10^{-37}$ @60 dps; box containment margin $1.5{\cdot}10^{-6}$ |
| H1: old breg harmless | `fixed/step_10_J192.py` | ${\rm breg\_abs} = 0.75 \ge B_{\rm abs} = 0.16666814$ (assert) |
| H2: lemmas + bound | `fixed/certify_H2.py` | $x^2|g|\le0.9999999995$; $2\mathrm{Re}(w)/|w|\ge1.4155$; AM-GM min $=0$; $I_0/{\rm bd}\approx0.27$ |
| H3: identity + boxes | `fixed/verify_mS.py` | resid $1.77{\cdot}10^{-31}$; upper-box margin $3.93{\cdot}10^{-21}$ (all $t$); lower-box margin $+1.31{\cdot}10^{-21}$ on $N\ge100$ |
| H3: lower box on domain | `fixed/certify_H3.py` | ${\rm Hp3}([0,0.13])\ge0.1538$ (arb); $H>0$ for $t<4.5447$ |
| sup $C$ $[100,400]$ | `fixed/step_18.py` | $194.1267251$ (29 cells) |
| sup $C$ sweep 1.5 | `fixed/planner_sweep_ratio.py 1.5` | $617118.8072$ (43 cells) |
| sup $C$ sweep 1.0002 | `fixed/planner_sweep_ratio.py 1.0002` | $254.2275176$ (85181 cells, 434 s) |
| final cell $C$ | `fixed/step_10_J192.py` | $47.0701250$ |
| theorem truth | `fixed/verify_theorem_truth.py`, D13r step_36 | $|M^{\rm main}-{\rm lead}|/{\rm sc} \le 1.07$ at $N=100,400,1891,10^4,10^5$ |

## Rerun commands

```bash
cd D13s_fix_recertify/D13s_work/fixed
python step_18.py                      # [100,400], 29 cells      -> sup C = 194.1267251   (~0.2 s)
python planner_sweep_ratio.py 1.5      # [400,1e10], 43 cells     -> sup C = 617118.8072   (~0.3 s)
python planner_sweep_ratio.py 1.0002   # [400,1e10], 85181 cells  -> sup C = 254.2275176   (~7.5 min)
python step_10_J192.py                 # final cell N >= 1e10     -> C     = 47.0701250    (~0.1 s)
```

Files: `fixed/step_38.py` (B2 + H2 fixes), `fixed/step_18.py`, `fixed/planner_sweep_ratio.py`,
`fixed/step_10_J192.py` (B1 + H1 + H2 + TLfN fixes), verifications `fixed/verify_A1.py`,
`fixed/verify_Bh.py`, `fixed/verify_zbox.py`, `fixed/verify_Btrue.py`, `fixed/certify_H2.py`,
`fixed/certify_H3.py`, `fixed/verify_TM_final.py`, `fixed/verify_TM_mono.py`,
`fixed/verify_theorem_truth.py`, `fixed/verify_mS.py`; buggy baselines preserved in `buggy/`.
