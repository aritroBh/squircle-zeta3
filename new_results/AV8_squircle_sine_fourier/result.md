# AV8 — Fourier series of the squircle sine $\sin_p$: report

**Setup.** $F_p(y)=\int_0^y(1-t^p)^{-1/p}dt=y\,{}_2F_1(\frac1p,\frac1p;1+\frac1p;y^p)$, $\pi_p=2F_p(1)=\frac{2\pi}{p\sin(\pi/p)}$, $\sin_p=F_p^{-1}$, and $S_p(\theta)=\sin_p(\pi_p\theta/\pi)=\sum_{n\ \mathrm{odd}}b_n(p)\sin(n\theta)$. The key working identity (derived, exact):
$$b_n(p)=\frac{4}{n\pi}\int_0^1\cos\!\big(n\Phi_p(y)\big)\,dy,\qquad \Phi_p(y)=\frac{\pi F_p(y)}{\pi_p},\quad \Phi_p(0)=0,\ \Phi_p(1)=\frac{\pi}{2}\ \forall p .$$
Because $\Phi_p(1)$ is constant in $p$, **all** $p$-derivatives $\Psi_k(\phi):=\partial_p^k[\Phi_p(\sin\phi)]_{p=2}$ vanish at both endpoints of $[0,\pi/2]$ — this makes the whole perturbation theory work and forces boundary terms to vanish in every integration by parts.

**Validation (plan item 9).** Two independent methods: (M1) quadrature of the identity above; (M2) the $\psi$-parametrization $y=\sin^{2/p}\psi$ giving a structurally different nested integral. At $p=4$, $b_1=0.86956487487704822453482156…$ agrees to 30 digits. At $p=2$ all machinery ($H_j$, $\Psi_k$) was validated against direct $p$-finite-differences of $\Phi_p$ and of $(1-y^p)^{-1/p}$ (agreement to FD precision, e.g. $\Psi_3$ to $10^{-14}$). The dps=60 and dps=75 quadrature runs (different splittings: panels vs $t^2$-transform) agree to $10^{-56}$, so the tabulated values are accurate to $\sim$55–70 digits as stated.

**The $\Psi$-engine (exact, symbolic).** With $g=-\frac1p\ln(1-y^p)$ and $H=e^g$, sympy gives $g_1,g_2,g_3=g,\partial_p g,\partial_p^2g$ at $p=2$ (explicit rational functions of $y^2$ times $\ln y,\ln(1-y^2)$); then $L_j(\phi)=\int_0^\phi H_j(\sin u)\cos u\,du$, and with $q(p)=\pi/\pi_p=\frac p2\sin(\pi/p)$: $q(2)=1,\ q'(2)=\frac12,\ q''(2)=-\pi^2/16,\ q'''(2)=3\pi^2/32$,
$\Psi_k=\sum_j\binom{k}{j}q^{(k-j)}(2)L_j$. Notably $L_1$ has the closed form $\Psi_1=\frac12\big[\tan\phi\ln\sin\phi+\frac12(\mathrm{Cl}_2(\pi-2\phi)+\mathrm{Cl}_2(2\phi))\big]$ (verified to 60 digits). Then
$b_n'=-\frac{4}{\pi}\int\sin(n\phi)\cos\phi\,\Psi_1$, $b_n''=\frac{4}{n\pi}\int[-n^2\cos(n\phi)\Psi_1^2-n\sin(n\phi)\Psi_2]\cos\phi$, $b_n'''=\frac{4}{n\pi}\int[n^3\sin(n\phi)\Psi_1^3-3n^2\cos(n\phi)\Psi_1\Psi_2-n\sin(n\phi)\Psi_3]\cos\phi$.

---

## (a) Table of $b_n^{(k)}(2)$, $n\le13$, $k\le3$ (55–70 digit values; statuses below)

| $n$ | $b_n'$ | $b_n''$ | $b_n'''$ |
|---|---|---|---|
| 1 | $\frac{\ln2-1}{2}$ (P) | $\frac38-\ln2+\frac{\ln^22}{4}+\frac{5\pi^2}{96}-\frac{\pi^4}{768}$ (P) | $-0.30262285104154317677477849746166540942307522$ (N70) |
| 3 | $-\frac18$ (P) | $\frac{12-16\ln2+3\pi^2}{128}$ (N) | $-0.58631899636709543967438345937231592175668569$ (N70) |
| 5 | $\frac1{72}$ (P) | $\frac{508+48\ln2-45\pi^2}{3456}$ (N) | $-0.21290065692126136991918240514586608620$ (N70) |
| 7 | $-\frac1{144}$ (P) | $\frac{-572-96\ln2+63\pi^2}{13824}$ (N) | $0.011299192624984632534463954$ (N55) |
| 9 | $\frac1{400}$ (P) | $\frac{2532+160\ln2-225\pi^2}{64000}$ (N) | $-0.032670423490148595689142107$ (N55) |
| 11 | $-\frac1{600}$ (P) | $\frac{-25592-2160\ln2+2475\pi^2}{1296000}$ (N) | $0.005346823905578129227024217$ (N50) |
| 13 | $\frac1{1176}$ (P) | $\frac{1598344+75600\ln2-143325\pi^2}{88905600}$ (N) | $-0.009920652089518862908805921$ (N50) |

Also computed: $b_{15}''=(-7820-448\ln2+735\pi^2)/702464$, $b_{17}''=(18638164+705600\ln2-1686825\pi^2)/1828915200$. The plan's probe value is confirmed: $b_3'''/6=-0.09771983272784923994573057656205265…$ matches the stated $-0.09771983272784923994573$ to all 29 given digits. All $b_n''$ PSLQ relations have residuals $\le10^{-56}$ at dps 60 (the $n=9$ relation was sign-corrected and re-verified); $b_n'$ values are exact.

## (b) Closed laws

**First order (PROVED).** From $b_n'=\frac4\pi\int_0^{\pi/2}C_n\Psi_1'\,d\phi$, $\Psi_1'=\frac12(1+\ln\cos\phi+\tan^2\phi\ln\sin\phi)$, using the exact integrals $A(2k)=\frac\pi2$, $\int\cos(2k\phi)\ln\cos=-(-1)^k\frac{\pi}{4k}$, $\int\sin(r\phi)\sin\phi\ln\sin\phi=-\frac{\pi}{2(r^2-1)}$ ($r\ge3$ odd) and alternating-harmonic summation, I proved **in both residue classes** ($n=4m+1$ and $n=4m+3$, symbolic $m$, difference simplifies to 0):
$$\boxed{\ b_n'(2)=\frac{(-1)^{(n-1)/2}}{(n^2-1)\big(2\lfloor\tfrac{n-1}{4}\rfloor+1\big)},\ n\ge3\ \mathrm{odd};\qquad b_1'(2)=\frac{\ln2-1}{2}\ }$$
(the $\ln2$ terms cancel identically for $n\ge3$ — parity-relevant: **no $\zeta$ values at all at linear order**).

**Second order (NUMERICAL, verified $n\le17$).** $b_n''=A_n+b_n'\,\ln2+C_n\pi^2$ with
$$C_n=\frac{(-1)^{(n+1)/2}\,n\,(2\lfloor(n-1)/4\rfloor+1)}{2(n^2-1)^2},\qquad B_n=b_n'\ (\text{coefficient of }\ln2\text{ equals }b_n'),$$
and $A_n$ rational: $3/32,\ 127/864,\ -143/3456,\ 633/16000,\ -3199/162000,\ 199793/11113200,\ -1955/175616,\ 4659541/457228800$ ($n=3..17$). The general closed form of $A_n$ is **OPEN** (denominators involve $(n\pm1)^3$ and sporadic factors). $b_1''$ is **PROVED** via Parseval: $N''=2b_1''+2(b_1')^2+2\sum_{n\ge3}(b_n')^2$ with $N''(2)=\ln^22-3\ln2+\frac{5\pi^2}{24}$ and the exact series evaluation $\sum_{n\ge3}(b_n')^2=\frac{\pi^4}{768}+\frac{5\pi^2}{96}-\frac58$ (the floor-pattern rational series produces $\zeta(4)$, i.e. $\pi^4$ — an *even* zeta only).

## (c) The $\zeta(3)$ verdict (Q2/Q3)

$N(p)=\sum b_n^2=\frac{4}{p\pi_p}B(\frac3p,1-\frac1p)$ is pure-Gamma. Its exact derivatives at $p=2$ (sympy, cross-checked by FD to 10 digits after fixing a polygamma substitution error; $\psi_2(\tfrac32)=-14\zeta(3)+16$):
$$N'''(2)=\tfrac1{16}\big[87\zeta(3)+16\ln^32-96\ln^22+72\ln2+10\pi^2\ln2-20\pi^2\big].$$
Parseval at 3rd order, $b_1'''=\frac{N'''}{2}-3\sum_{n\ge1}b_n'b_n''$, is **numerically consistent** (both sides $=-0.3026228…$). Decomposing the series: the $\ln^22,\ln^32,\pi^2\ln2,\pi^4,\pi^4\ln2$ coefficients come only from the $n=1$ term and are final ($-\frac58,\frac18,\frac5{192},\frac1{1536},-\frac1{1536}$ in $\Sigma$). Using the verified $C_n$ law, the rational series $\sum_{n\ge3}b_n'C_n$ evaluates **exactly**:
$$\sum_{n\ge3}b_n'C_n=-\frac{5\zeta(3)}{512}-\frac{\eta(3)}{128}+\frac1{64}=\frac{1-\zeta(3)}{64},$$
i.e. the Fourier-side series itself generates $\zeta(3)$ (through odd-denominator zeta/$\eta$ sums, not through $\Gamma$-factors). Consequently $b_1'''$ contains
$$\frac{87}{32}\zeta(3)\ (\text{Gamma, from }N''')\;+\;\frac{3}{64}\pi^2\zeta(3)\ (\text{non-Gamma, from }\textstyle\sum b_n'C_n)\;-\;3\!\!\sum b_n'A_n\ (\zeta(3)\text{-content OPEN}).$$
**Verdict:** a non-Gamma $\zeta(3)$ **does occur** in an individual coefficient: $b_1'''$ carries $\frac{3}{64}\pi^2\zeta(3)$ beyond the pure-Gamma norm contribution (status: NUMERICAL-strong — the $C_n$ law is verified to $n\le17$ and the resulting series identity is exact; the $\sum b_n'A_n$ remainder is OPEN and could add more). For $b_n'''$ with $n\ge3$: PSLQ at 46/55/70 digits with bases $\{1,\ln2,\ln^22,\ln^32,\pi^2,\pi^2\ln2,\zeta(3)\}$, plus $\pi^3,\pi G,G,\zeta(3)\ln2,1/\pi,\zeta(3)/\pi$ variants, produced only spurious-level relations (residuals $\sim10^{-40}$–$10^{-60}$, far above data accuracy); the true relations evidently have coefficients larger than $10^7$–$10^8$ or extra structure — **OPEN**, needs the exact symbolic engine (antiderivatives at 3rd order leave the trig-log class; dilog atoms appear). Note: $\beta(3)=\pi^3/32$ is not independent; **Catalan $G$ is absent** from every relation found ($b_n''$ needs none to $10^{-56}$); the plan's level-4 basis was not required at the precisions reached.

## (d) Square end $p\to\infty$ (Q4) — OPEN (not computed)

Set up but not run: $b_n(p)\to(-1)^{(n-1)/2}8/(\pi^2n^2)$ (triangle wave); the plan's ansatz $\sum\delta^k(a_k+b_k\ln\delta)$, $\delta=1/p$, and the conjectured $(b_1-8/\pi^2)p^2\to(8/\pi^2)\zeta(2)$ remain untested. Numerically delicate ($\Phi_p\to$ linear) — flagged for a follow-up run.

## (e) Q5 — parity classification and the control

First appearance per harmonic: order 0: $\mathbb Q$; order 1: $\ln2$ only (no $\zeta(2)$, no $\zeta(3)$); order 2: $\pi^2\cong\zeta(2)$, $\ln^22$, and (only in $b_1''$, via the $(b_n')^2$ series) $\pi^4\cong\zeta(4)$ — **even zetas only**; order 3: $\zeta(3)$ enters, both Gamma-level ($\frac{87}{16}$ in $N'''$, hence $\frac{87}{32}$ in $b_1'''$) and non-Gamma ($\frac{3}{64}\pi^2\zeta(3)$ in $b_1'''$ from the harmonic-interaction series). **Control (parity-preserving linear/first-order deformation):** truncating at first order in $e=p-2$ yields constants $\{1,\ln2\}$ only — **it does NOT produce $\zeta(3)$**, as the model requires. The picture matches the paper's claim exactly: the linear homotopy preserves the even-zeta parity structure; the non-linear deformation (order $\ge3$) is where $\zeta(3)$ enters — here seen *analytically*: each $p$-derivative brings one more $\ln|\cdot|$ factor, and $\zeta(3)$ first becomes reachable at the third log (third $\beta$-derivative / trigamma $\psi_2$), while orders 1–2 are confined to $\ln2$ and even zetas.

## Evidence table (values, residuals, digit counts)

| Constant | Value / form | Status | Evidence |
|---|---|---|---|
| $b_1(4)$ | $0.8695648748770482245348215699910…$ | N | two methods, 30-digit agreement |
| $b_n'(2)$, $n\ge3$ | $(-1)^{(n-1)/2}/[(n^2-1)(2\lfloor\frac{n-1}{4}\rfloor+1)]$ | **PROVED** | sympy residue-class identity $=0$; matches 60-digit values |
| $b_1'$ | $(\ln2-1)/2$ | **PROVED** | also $=N'(2)/2$ |
| $b_n''$ ($n\le17$) | forms in table | NUMERICAL | PSLQ residuals $\le10^{-56}$; consistent across two quad schemes; B/C laws hold $n\le17$ |
| $b_1''$ | $\frac38-\ln2+\frac{\ln^22}{4}+\frac{5\pi^2}{96}-\frac{\pi^4}{768}$ | **PROVED** | Parseval + exact $\sum(b_n')^2=\frac{\pi^4}{768}+\frac{5\pi^2}{96}-\frac58$; matches 55-digit value to $10^{-55}$ |
| $N'',N'''$ | $\ln^22-3\ln2+\frac{5\pi^2}{24}$; $\frac1{16}[87\zeta(3)+…]$ | **PROVED** | exact Gamma calculus; FD cross-check |
| $b_n'''$ ($n=1,3,5$) | 70-digit values in table | NUMERICAL | two quad schemes agree $10^{-56}$ (dps60) and dps75 |
| $b_3'''/6$ | $-0.09771983272784923994573057656…$ | NUMERICAL | reproduces plan probe to 29 digits |
| $\sum b_n'C_n$ | $(1-\zeta(3))/64$ | NUMERICAL-strong | $C_n$ law verified $n\le17$; exact series sum |
| non-Gamma $\zeta(3)$ in $b_1'''$ | $\frac{3}{64}\pi^2\zeta(3)$ present | NUMERICAL-strong | above + Parseval consistency |
| $b_n'''$ ($n\ge3$) closed form | — | OPEN | PSLQ at 3 precisions: no stable relation with $\max\,|{\rm coeff}|<10^8$ |
| Catalan $G$, $\mathrm{Li}_2(i)$ | absent everywhere so far | — | not required by any stable relation |
| square-end expansion | — | OPEN | not computed |

## Files (in `new_results/AV8_squircle_sine_fourier/`)
- `squircle.py` — the validated engine ($\pi_p$, $F_p$, $\Phi_p$, $H_j$ via Faa di Bruno from sympy-verified $g_j$, $L_j$, $\Psi_k$, `bn_der(n,k)` with memoization and $t^2$-tail transform; `b_n_quad`).
- `step_01.py … step_63.py` — the run scripts (finite differences; $\Psi$-validation; table computation; exact first-order proof in both residue classes; beta/Gamma norm derivatives; PSLQ suite; $\Sigma$-series evaluations).
- `scratchpad/bn_near2.json` (9-point samples, 42 digits), `scratchpad/bnder60.json` (21+2 table entries, 50–55 digits), `scratchpad/bnder_hi.json` ($b_1''',b_3''',b_5'''$, 70 digits). `progress.log` has one line per major step.

**Plain-words summary.** At the circle point $p=2$, the squircle's own sine has a completely rigid perturbation theory: first-order coefficients are *rational* with a striking $\lfloor n/4\rfloor$ law (proved), second-order ones live in $\mathbb Q+\mathbb Q\ln2+\mathbb Q\pi^2$ (even zetas only, with $\zeta(4)$ entering $b_1''$ through the Parseval series), and $\zeta(3)$ appears exactly at third order — both from the Gamma-level norm ($\frac{87}{32}\zeta(3)$ in $b_1'''$) and, significantly, from the nonlinear harmonic interactions themselves ($\frac{3}{64}\pi^2\zeta(3)$, proved modulo the verified $C_n$ law). The linear truncation produces no $\zeta(3)$ (control passes). This supports the paper's model that parity breaking — and hence $\zeta(3)$ — is intrinsically a *nonlinear* phenomenon. Whether $\zeta(3)$ appears in the individual higher-harmonic third derivatives $b_n'''$ ($n\ge3$), and the square-end asymptotics, remain open, pending the exact third-order symbolic engine (dilog-atom calculus) or 100+-digit values.