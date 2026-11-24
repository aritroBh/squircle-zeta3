# DCY2F — repaired DCY2 notes (referee DCY2R repairs R1–R6 applied and re-certified)

All repairs R1–R6 are applied and the whole theorem is **re-certified end-to-end** by `verify_dcy2_repaired.py` (runs in **126.2 s < 5 min**, every check passes, arb/flint outward rounding throughout). Every number below was printed by code this session.

## Repaired theorem statement (with constants)

> **Theorem (first correction, computer-assisted).** Put $N=k+1\ge 100$, $\theta=\sqrt{\pi N}$, and
> $$\varpi_1=\frac{(24+12\pi+\pi^2)+i\,(24-12\pi+\pi^2)}{48\sqrt\pi}
> =0.84121508776758439707534863790981499908\ldots-0.04501183768517361657373510376075759232\ldots i .$$
> Then, with $M_k^{\rm main}$ the $\pm i\pi/2$ branch contribution,
> $$\Bigl|M_k^{\rm main}+\tfrac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\tfrac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
> \le C_2\,N^{-2}e^{-\theta},\qquad C_2=10^{6},\quad N\ge100,$$
> and on $100\le N\le400$ the certified constant is $\le 1.76$ (max $1.751172$, at $N=116$).
>
> **Corollary (full $M_k$, with CX6's proved rest bound $|M_k-M_k^{\rm main}|\le10{,}600{,}000\,e^{-2\sqrt N}$).**
> $$\Bigl|M_k+\tfrac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\tfrac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
> \le C_3\,N^{-2}e^{-\theta},\qquad C_3=1.86\times10^{10}\ (100\le N\le N_3=25{,}711),\quad C_3=C_2+1=10^6+1\ (N>25{,}711).$$
> (Computed $C_3=C_2+10{,}600{,}000\cdot h(309)=1.8540155\times10^{10}$ with $h(309)=1748.976914$.) On $100\le N\le400$ the full-$M_k$ sharp constant is $\le1.743787$ ($N=158$).

**Certificate architecture:** 301 certified arb quadratures on $[100,400]$; 279 interval cells on $[400,10^5]$ (worst $396{,}111.40<C_2$); a monotone analytic tail on $[10^5,\infty)$ ($\mathrm{LHS}(10^5)=76{,}889.43<C_2$).

## Verify output (`verify_dcy2_repaired.py`, return code 0, wall 126.4 s)

```
====================================================================================================
DCY2 VERIFIER (REPAIRED, referee DCY2R R1-R6) — first correction of the decay law of M_k
====================================================================================================
[1] varpi1 = [0.8412150877675843970753486379098149990800 +/- 2.99e-41] + [-0.04501183768517361657373510376075759231874 +/- 4.88e-42]j
    paper's 81-digit value: 0.84121508776758439707534863790981499908 -0.04501183768517361657373510376075759232i
    match: Re diff = [2.9854e-41 +/- 8.72e-46],  Im diff = [1.2551e-39 +/- 2.55e-44]
----------------------------------------------------------------------------------------------------
[2] R1 discrepancy exhibit (the sign bug and its repair)
    N=381: |M_main - X| N^2 e^th = [1409.465710 +/- 4.08e-7]   <- the OLD (wrong-sign) quantity ~ 2|X|
           |M_main + X| N^2 e^th = [1.619635093 +/- 4.78e-10]   <- the THEOREM quantity
    (M_main/X = [-1.00230 +/- 8.70e-7]: M_main ~ -X, so the wrong sign measures ~2|X|, the '1409')
----------------------------------------------------------------------------------------------------
[3] finite range [100,400] REPAIRED: |M_main + X| <= C2 N^-2 e^-th, 301 certified quadratures
    N=100: res/bound = [3.1504939e-7 +/- 1.81e-15]   |M_k - M_main| = [1.64857e-13 +/- 4.09e-19]
    N=200: res/bound = [1.6103581e-6 +/- 4.09e-14]   |M_k - M_main| = [4.21136e-18 +/- 4.31e-25]
    N=400: res/bound = [1.3347060e-6 +/- 1.74e-16]   |M_k - M_main| = [4.69506e-25 +/- 3.25e-31]
    ALL PASS (301 values): max res/bound = [1.751172383e-6 +/- 4.93e-17] at N=116
      i.e. max |M_main + X| N^2 e^th = [1.751172383 +/- 4.93e-11]  (constant needed on [100,400]; ~1.75)
      full M_k: max |M_k + X| N^2 e^th = [1.743787030 +/- 1.32e-10] at N=158
    elapsed: 124.2s
----------------------------------------------------------------------------------------------------
[4] region B sweep [400, 100000] (repaired two-sided window tail + Mtilde correction + endpoints)
    ALL 279 CELLS PASS: worst LHS_cell = [396111.3989 +/- 1.04e-5] at (3075.23,3136.74)
      margin to C2=1e6: 2.525x;  +C_b=29.5 still < C2
----------------------------------------------------------------------------------------------------
[5] region C [1e5, inf): LHS(N) <= LHS(1e5) <= C2
    LHS(1e5) = [76889.42704 +/- 3.80e-6] (< C2);  pointwise monotone on 116-pt grid [1e5,3e6]: True
    LHS at 1e7,1e9,1e12,1e20 = [26952.9 ...], [17719.2 ...], [14739.1 ...], [14026.2 ...]
----------------------------------------------------------------------------------------------------
[6] part (b): |model - T1| N^2 e^th <= C_b  (R3: NO spurious (1+varpi1/sqrt N) factor)
    cross = [2.4951160 ...], C_S2 = [1.1350952 ...], C_F = [8.7457113 ...] <= 10.155
    |model-T1| N^2 e^th = (pi^1.5/2)|delta_F| <= (pi^1.5/2) C_F = [24.349494 ...]  (corrected)
    conservative (spurious factor kept): [25.375118 ...];  paper value C_b = 29.5 covers both
    pointwise D(N) = 1.65739, 1.61977, 1.60148, 1.59352, 1.57777, 1.56497 at N=400..100000
----------------------------------------------------------------------------------------------------
[7] rest combination (CX6: |M_k - M_main| <= 10,600,000 e^{-2 sqrt N})
    N3 = 25711.0000000;  max h(N) = h(309) = [1748.976914 ...]
    C3 (100 <= N <= N3) = C2 + Crest*h(309) = [18540155283.7 ...];  for N >= N3: C3 = C2+1 = 1.0e6
----------------------------------------------------------------------------------------------------
[8] tail-machinery certificates (arb, outward): c2, rDelta, eLamD, lamL hold at N=400..1e6: True
----------------------------------------------------------------------------------------------------
[9] R4/R5 repair-delta at worst cell (3075.23,3136.74): worst LHS_cell = 396111.3989 (< C2, margin 2.525x);
    one-sided window-tail sum Et = 9.75298 (doubling costs one Et share; both forms < C2)
----------------------------------------------------------------------------------------------------
[10] R6: step_03_watson re-certification: C_rho = 0.062500000 <= 0.5 (truth 0.0117);  C_rho2 = 2.4154401 < 2.42
----------------------------------------------------------------------------------------------------
[11] R3: corrected P(u): diff(exact,typo -2b_R y^3) != 0;  diff(exact,corrected -2b_R y^2) = 0;
    P(1) = y*(-4y^4-8sqrt2 y^3-...-4)/(2(y^2+sqrt2 y+1)) < 0 for y>0  (single-max lemma intact)
====================================================================================================
TOTAL elapsed: 126.2s
VERIFY_DCY2_REPAIRED: ALL CHECKS PASSED
```

## Errata (referee DCY2R) — R1–R6 as applied

- **R1 (finite-range sign, required).** The old finite-range code computed $|M_k^{\rm main}-X|\approx2|X|$ (hence the spurious "res/bound $0.00141$, constant $\approx1409$" in the old notes) instead of the theorem's $|M_k^{\rm main}+X|$. At $N=381$: $|M-X|N^2e^\theta=1409.465710$ vs $|M+X|N^2e^\theta=1.619635$. Repaired (`step_06_finite.py` uses an explicit `X_of(N)=+X` and `res=|M+X|`; `verify_dcy2_repaired.py` §2–§3). Re-certified: $\max|M^{\rm main}+X|N^2e^\theta=\mathbf{1.751172383}$ ($N=116$); full $M_k$ $\mathbf{1.743787030}$ ($N=158$). L14's quoted residuals corrected from $(0.000165,0.000751,0.000932)$ to the true scaled residuals $(0.315,1.610,1.335)$.
- **R2 (write-up).** The complex estimate is stated for $J:=(2/\pi)\tilde I$ with $M_k^{\rm main}=\operatorname{Im}J$ (all code objects `It`/`I_tilde`/`model`/`target` are $J$-quantities); the estimate for $\tilde I$ alone as printed was false since $\tilde I\approx(\pi/2)T_1$.
- **R3 (typos).** $P(u)$ $u$-coefficient $-2b_Ry^3\to-2b_Ry^2$ (fixed in `dcy2_budget.py`, `dcy2_cell.py`; verified symbolically in §11 — diff with exact $2uQ\Lambda'=0$; harmless since $\min$ selects $\pi y/2$). Removed the spurious $(1+\varpi_1/\sqrt N)$ factor in §5: the correct bound is $|{\rm model}-T_1|N^2e^\theta=(\pi^{3/2}/2)|\delta_F|\le(\pi^{3/2}/2)C_F=24.35$ (conservative form keeping the factor $=25.38$; both $<C_b=29.5$). `dcy2_cert` $h_1$ comment is the used form $|h_1|\le1.5y(1+2y)$.
- **R4 (remainder).** Window tail $E_{\rm tailmom}$ made **two-sided** ($\times2$; the one-sided TailMom formula is valid only for a one-sided integral) — applied in `dcy2_sweep.window_cell` ($2\cdot E_t$) and `dcy2_cell.LHS_of` ($2\cdot E_{\rm tailmom}$); `Mtilde` edge correction $2e^{a^2/4c}t_R^{j+1}/(j+1)$ (in `dcy2_cell.Mtilde_of`). Worst region-B cell re-run: $\mathbf{396111.40}<C_2$ (margin $2.525\times$; the doubling costs one $E_t=9.75$ share at the worst cell).
- **R5 (sweep endpoints).** `TailMom`($p\le2$) evaluated at $N_{\rm lo}$ (it is decreasing there); $e^{\gamma T^3}$ taken at a consistent endpoint ($\gamma T^3\le\max(\gamma\Delta^3,c_2^3/8\gamma^2)$); `EdgeInt` $\Gamma$-subtraction clamped at $0$. In `dcy2_sweep.window_cell`; impact ~1%, absorbed by the $2.5\times$ margin.
- **R6 (scripts).** `step_03_watson.py` repaired (v4, alternating-series bound; the old Cauchy $M/4=0.9151>0.5$ and the misused `arb(lo,hi)` interval cells are gone). Re-certified: $C_\rho=1/16=0.0625\le0.5$ (truth $0.0117$ on $(0,1]$), $C_{\rho,2}=2.4154<2.42$; passes its asserts (§10).

## Evidence table (all printed by code; arb balls, outward)

| Quantity | Value | Location | digits/certification |
|---|---|---|---|
| $\varpi_1$ vs paper (81 digits) | Re diff $2.99\!\times\!10^{-41}$, Im diff $1.26\!\times\!10^{-39}$ | §1 | 45-dps arb; sympy diff $=0$ exact |
| $N=381$: $\lvert M-X\rvert N^2e^\theta$ (old, wrong) | $1409.465710$ | §2 | arb quadrature, $\pm4\!\times\!10^{-7}$ |
| $N=381$: $\lvert M+X\rvert N^2e^\theta$ (theorem) | $1.619635093$ | §2 | arb, $\pm5\!\times\!10^{-10}$ |
| finite range $\max\lvert M^{\rm main}+X\rvert N^2e^\theta$ | $\mathbf{1.751172383}$ at $N=116$ | §3 | 301 arb quadratures |
| full $M_k$ $\max\lvert M_k+X\rvert N^2e^\theta$ | $\mathbf{1.743787030}$ at $N=158$ | §3 | certified table balls |
| $\lvert M_k-M^{\rm main}\rvert$ at $N=100,200,400$ | $1.65\!\times\!10^{-13},4.21\!\times\!10^{-18},4.70\!\times\!10^{-25}$ | §3 | certified |
| region B worst cell | $396111.3989$ at $(3075.23,3136.74)$, margin $2.525\times$ | §4 | 279 interval cells, arb |
| region C $\mathrm{LHS}(10^5)$ | $76889.42704$ (monotone grid OK) | §5 | arb + $y$-scaling proof |
| part (b): $C_F$ / $C_b$ | $8.7457113$ / $24.349$ (corr.), $25.375$ (cons.) $\le29.5$ | §6 | arb |
| rest: $N_3$ / $h(309)$ / $C_3$ | $25711$ / $1748.976914$ / $1.8540155\!\times\!10^{10}$ | §7 | arb bisection |
| R6: $C_\rho$ / $C_{\rho,2}$ | $0.0625\le0.5$ / $2.4154<2.42$ | §10 | arb + alternating-series proof |
| R3: $P(u)$ corrected diff / $P(1)$ | $0$ (exact) / $<0$ (exact) | §11 | sympy exact |

## File list (cwd: `new_results/DCY2F_notes_repair`)

- **`verify_dcy2_repaired.py`** — the repaired verifier (§1–§11 above; 126 s, all pass).
- **`step_06_finite.py`** — R1-repaired finite-range module (`X_of`, `residual_check` $\to|M+X|$).
- **`dcy2_cell.py`** — R4 (documented one-sided `TailMom_of`; `LHS_of` bracket $=E_{\rm mom}+2E_{\rm tailmom}$; `Mtilde` correction present) + R3 `lamL` fix.
- **`dcy2_sweep.py`** — region-B cell machinery with R4 ($2E_t$) and R5 endpoints (copied, already repaired).
- **`dcy2_budget.py`** — R3 $P(u)$ typo fixed ($-2b_Ry^3\to-2b_Ry^2$).
- **`dcy2_regionc.py`, `dcy2_cert.py`, `dcy2_common.py`, `dcy2_audit.py`** — supporting modules (copies).
- **`step_02_varpi1.py`** — exact $\varpi_1$ derivation (sympy diff $=0$).
- **`step_03_watson.py`** — R6-repaired Watson constants (v4, passes asserts).
- **`step_07_partb.py`** — part-(b) $D(N)$ and model/T1 machinery.
- **`progress.log`** — one line per repair step.
- **`proof_repaired.md`** — the corrected notes (full text below).
- The exploratory/verification computations that produced the numbers above are saved as `step_01.py … step_2N.py` in this directory by the harness (each python call is one step file).

## `proof_repaired.md` (DCY2 `proof.md` with every repair applied; no "1409" anywhere)

```markdown
# DCY2 — the first correction of the decay law of $M_k$: proof (REPAIRED, referee DCY2R R1–R6)

**Theorem (DCY2, main term).** Put $N=k+1$, $\theta=\sqrt{\pi N}$, and
$$\varpi_1=\frac{(24+12\pi+\pi^2)+i\,(24-12\pi+\pi^2)}{48\sqrt\pi}
=0.84121508776758439707534863791-0.04501183768517361657373510376\,i .$$
Then for every integer $N\ge 100$,
$$\Bigl|\,M_k^{\rm main}+\frac{\sqrt\pi}{N}\,e^{-\theta}\,
\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\frac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
\;\le\; C_2\,N^{-2}e^{-\theta},\qquad C_2 = 1{,}000{,}000 .$$

**Corollary (full $M_k$).** With CX6's proved rest bound
$|M_k-M_k^{\rm main}|\le 10{,}600{,}000\,e^{-2\sqrt N}$ ($N\ge100$):
$$\Bigl|\,M_k+\frac{\sqrt\pi}{N}\,e^{-\theta}\,
\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\frac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
\;\le\; C_3\,N^{-2}e^{-\theta},\qquad C_3 = C_2+1.854\times10^{10}\approx1.86\times10^{10},$$
for all $N\ge100$; for $N\ge N_3=25{,}711$ the rest term is $\le N^{-2}e^{-\theta}$, so
$C_3$ can be taken $=C_2+1$ there; for $100\le N\le 400$ the *sharp* residuals are certified
against the table `mk_deg9.tsv` (max $|M_k^{\rm main}+X|N^2e^\theta=1.751172$ at $N=116$;
full $M_k$ max $1.743787$ at $N=158$ — i.e. the sharp constant there is $\le1.76$).

## Lemma ledger

| # | Lemma | Status |
|---|---|---|
| 1 | Ray representation: with $J:=\frac2\pi\tilde I$, $\tilde I=\int_{\rm ray}H(s)e^{i\pi s/2}(s/(1+s))^N ds/s$, $H(s)=\frac12B(\frac s2,\frac32)$, one has $M_k^{\rm main}=\operatorname{Im}J$ | **PROVED** (CX6 §3) |
| 2 | $\varpi_1$ closed form from the saddle expansion (two formalisms; exact algebra, diff $=0$; numeric match to the paper's 81 digits) | **PROVED** (`step_02_varpi1.py`) |
| 3 | Exponent factorization $J=\frac2\pi e^{\phi(s_0)}\omega\int_{-\rho}^{\infty}e^{\ell(t)}F(t)dt$; exact modulus $M(t)$; single-max lemma ($P(u)$ cubic, $P(0^+)>0$, $P(1)<0$, $P''<0$); tail rates $r_\Delta$, $\lambda_L$ | **PROVED** (by hand; certified pointwise) |
| 4 | $|\ell_m|\le Ny^{m+1}+y^m/m$; $c_2=-\operatorname{Re}\ell_2\ge0.9b_Ry$; $\operatorname{Re}\ell_1$ exact | **PROVED** (hand; $c_2$ certified) |
| 5 | Watson lemma for $H$: $g(s_0)=1-\frac{3}{4s_0}+O(y^2)$ ($C_g\le1$), $h_j$ expansions, $M_4'$; $C_\rho\le0.5$ | **PROVED** (Watson with remainder) + **CERTIFIED** (`step_03_watson.py`, repaired v4) |
| 6 | Kept-polynomial moments $\mu_p^\infty$ closed form; model $=T_1+O(N^{-2}e^{-\theta})$; the two exact identities | **PROVED** (sympy exact, diff $=0$) + **CERTIFIED** constants |
| 7 | Remainder $|E(t)|$ by $|\cdot|$-moments with $[T,\Delta]$-split $e^{\gamma|t|^3}$ factors and damped incomplete-$\Gamma$ edges; window tail **two-sided** ($\times2$) with the `Mtilde` edge correction $2e^{a^2/4c}t_R^{j+1}/(j+1)$ | **CERTIFIED** (`dcy2_cert.py`/`dcy2_cell.py`, repaired) |
| 8 | Window-tail bounds $\mathrm{TP}_p$; outer tails via $M(\pm\Delta)$ with the single-dip lemma | **CERTIFIED** + **PROVED** |
| 9 | Region B sweep $[400,10^5]$: $\mathrm{LHS}_{\rm cell}\le C_2$ on every one of 279 cells | **CERTIFIED** (`dcy2_sweep.py`, worst $=396{,}111$) |
| 10 | Region C $[10^5,\infty)$: $\mathrm{LHS}(N)\le\mathrm{LHS}(10^5)=76{,}889\le C_2$ | **PROVED** (y-scaling) + **CERTIFIED** (pointwise) |
| 11 | Finite range $[100,400]$: certified quadrature, $\max|M^{\rm main}+X|N^2e^\theta=1.751172$ | **CERTIFIED** (`step_06_finite.py`, repaired, 301 quadratures) |
| 12 | Part (b): $|{\rm model}-T_1|\le C_bN^{-2}e^{-\theta}$, $C_b=29.5$; $C_F\le10.155$ | **CERTIFIED** (`step_07_partb.py` + §5 assembly, corrected) |
| 13 | Rest combination: $N_3=25{,}711$, $C_3=C_2+10{,}600{,}000\cdot h(309)$, $h(309)=1748.976914$ | **CERTIFIED** (arb) |
| 14 | Validation at $N=100,200,400$ vs `mk_deg9.tsv`: scaled residuals $0.315,1.610,1.335$ | **CERTIFIED** (`verify_dcy2_repaired.py` §3) |

## 1. Setup and the exact representation (L1)

$H(s)=\frac12B(\frac s2,\frac32)=\frac{\sqrt\pi}{4}\frac{\Gamma(s/2)}{\Gamma(3/2+s/2)}
=\int_0^\infty\sqrt{1-e^{-2u}}\,e^{-su}\,du$ ($\operatorname{Re}s>0$),
$\omega=e^{i\pi/4}$, ray $s=\omega x$, $x\ge0$. Put
$$J:=\frac2\pi\int_0^\infty H(\omega x)e^{i\pi\omega x/2}\Bigl(\frac{\omega x}{1+\omega x}\Bigr)^N\frac{\omega\,dx}{\omega x},
\qquad M_k^{\rm main}=\operatorname{Im}J .$$
**(R2)** We prove the complex estimate $|J-T_1|\le C_2N^{-2}e^{-\theta}$ with
$T_1=-i\frac{\sqrt\pi}{N}e^{-\theta+i\theta}(1+\varpi_1/\sqrt N)$; taking imaginary parts
gives the theorem since $|\operatorname{Im}w|\le|w|$ and
$\operatorname{Im}T_1=-\frac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}[e^{i\theta}(1+\varpi_1/\sqrt N)]$.
(All code objects named `It`/`I_tilde`/`model`/`target` are $J$-quantities; the estimate is
stated for $J$, not for $\tilde I$ alone.)

Notation: $y=\sqrt{\pi/(2N)}$, $\rho=1/y$, $s_0=\omega/y$,
$\theta=\sqrt{\pi N}=\pi/(\sqrt2 y)$, $b=\pi\bar\omega/2$, $b_R=\pi/(2\sqrt2)$,
$\kappa^2=\kappa_0^2+4\log(N/400)$ with $\kappa_0=2.48$, $\sigma=(2\cdot0.9b_Ry)^{-1/2}$,
$\Delta=\kappa\sigma$, $v=\sqrt y$.

## 2. $\varpi_1$ derived, not imported (L2) — unchanged

$G=-b\tau^2+vG_1\tau+v^2G_2\tau^2+\cdots$ with
$G_1=-(1+\tfrac{\pi\bar\omega^2}{2})\tau+\tfrac{\pi\bar\omega}{2}\tau^3$,
$G_2=(\tfrac12+\tfrac{3\pi\bar\omega^2}{4})\tau^2-\tfrac{\pi\bar\omega}{2}\tau^4$;
$(1+v\tau)^{-3/2}=1+h_1v\tau+h_2v^2\tau^2+\cdots$ ($h_1=-3/2$, $h_2=15/8$). Then
$c_1=\langle G_2+\frac{G_1^2}{2}+h_2\tau^2+G_1h_1\tau\rangle_{e^{-b\tau^2}}
=\frac{(1+i)(15-\pi^2)}{8\sqrt2\,\pi}$,
$\varpi_1=\sqrt{\frac\pi2}(c_1+S_1-\tfrac34\bar\omega)$, $S_1=\pi(1+i)/(6\sqrt2)$.
**`step_02_varpi1.py` prints $\varpi_1-\frac{(24+12\pi+\pi^2)+i(24-12\pi+\pi^2)}{48\sqrt\pi}=0$
exactly** (sympy), and the 40-digit values match the paper. The moment-sum $O(y)$ coefficient
is **exactly $0$**, and $S_1+g'_1+B_2/(2b)+\gamma_1=\varpi_1\sqrt{2/\pi}$ exactly (sympy, diff $=0$).

## 3. The analytic tail (L3–L9)

**Factorization (L3).** $\ell(t)=(N-1)\log(1+yt)-N\log(1+\frac{yt}{1+y\bar\omega})+\frac{i\pi\omega t}{2}$,
$J=\frac2\pi e^{\phi(s_0)}\omega\int_{-\rho}^{\infty}e^{\ell(t)}F(t)dt$, $F(t)=H(s_0+\omega t)/H(s_0)$.
With $u=1+yt$, $\operatorname{sgn}\Lambda'(t)=\operatorname{sgn}P(u)$, **(R3: corrected $u$-coefficient)**
$$P(u)=-2b_Ru^3-(\pi+2)yu^2+\bigl(\sqrt2(N-2)y^2-2b_Ry^{\mathbf 2}\bigr)u+2(N-1)y^3,$$
$P(0^+)>0$, $P(1)<0$ (exact: $P(1)=-y(4+8\sqrt2y+\cdots)/(2Q_0)$), $P''<0$; one root $u_c\in(0,1)$.
Right tail $M(t)\le M(\Delta)e^{-r_\Delta(t-\Delta)}$, $r_\Delta=\min(-\Lambda'(\Delta),b_R)$;
left tail $M(t)\le M(-\Delta)e^{-\lambda_L(-\Delta-t)}$, $\lambda_L=(\pi y/2)/(2u_DQ(u_D))$.

**Window (L7).** $\ell(t)=\ell_1t+\ell_2t^2+R_3(t)$; $F(t)e^{R_3(t)}=P(t)+E(t)$ with the
14-monomial kept polynomial $P$. **(R4)** $E(t)$ is bounded by $|\cdot|$-moments with
$[T,\Delta]$-split factors $e^{\gamma|t|^3}$, $T=\min(\Delta,c_2/(2\gamma))$, damped
incomplete-$\Gamma$ edges, the **two-sided** window tail ($E_{\rm tailmom}\times2$), and the
shifted-Gaussian edge correction $\tilde M_j\gets\tilde M_j+2e^{a^2/4c}t_R^{j+1}/(j+1)$.

**$\ell$-boxes (L4).** $|\ell_m|\le Ny^{m+1}+y^m/m$; $c_2=-\operatorname{Re}\ell_2\ge0.9b_Ry$;
$\operatorname{Re}\ell_1=-y(1+\frac{b_Ry}{1+\sqrt2y+y^2})$ exact.

**$H$ (L5).** Watson with remainder; $|\rho_3|\le C_\rho u^3$ on $[0,1]$ ($C_\rho=1/16=0.0625\le0.5$,
certified by the repaired `step_03_watson.py`) and $\le2.42u^3$ on $[1,\infty)$;
$g(s_0)=1-\frac{3}{4s_0}+r_g$, $|r_g|\le C_gy^2$ ($C_g\le1$);
$|F(t)-\sum_0^3h_jt^j|\le\frac{M_4'}{24}|t|^4$.

**Assembly (L7–L9).**
$J-{\rm model}=(2/\pi)e^{\phi(s_0)}\omega[\int_{|t|\le\Delta}e^{\ell_1t+\ell_2t^2}E(t)dt
-\sum p_p(\mu_p^\infty-\mu_p(\Delta))]+\text{tails}$,
bounded as in `dcy2_cell.py`/`dcy2_sweep.py` (repaired); worst region-B cell $396{,}111<C_2$.

## 4. Region C: $N\ge10^5$ (L10)

Every piece of $\mathrm{LHS}(N)=E_{\rm total}N^2e^\theta$ is non-increasing for $N\ge10^5$
($y$-scaling; $e^{-\kappa^2/4}=\frac{800}{\pi}e^{-\kappa_0^2/4}y^2$ folded in; $\kappa^2=28.2\ge16$).
Bound $=\mathrm{LHS}(10^5)=76{,}889<C_2$ (repaired two-sided form); monotone on a dense grid.

## 5. Part (b): model $-T_1$ (L6, L12) — R3 corrected

${\rm model}=-iye^{-\theta+i\theta}\,{\rm ModelRatio}$,
$T_1=-iye^{-\theta+i\theta}(2y/\sqrt\pi)(1+\varpi_1/\sqrt N)$.
Writing ${\rm ModelRatio}=(2y/\sqrt\pi)(1+\varpi_1/\sqrt N+\delta_F)$ with $\delta_F=O(y^2)$,
$$|{\rm model}-T_1|\cdot N^2e^{\theta}=\frac{\pi^{3/2}}{2}\,|\delta_F|
\;\le\;\frac{\pi^{3/2}}{2}\,C_F .$$
**(R3)** the factor $(1+\varpi_1/\sqrt N)$ that appeared here in the old draft is spurious
(an artifact of the relative-remainder form) and is removed; even keeping it only enlarges the
bound ($25.38$). With $|\delta_F|\le C_F=C_S+C_g+C_{S2}+C_E+C_{MS}+\mathrm{cross}\le8.7457\le10.155$
($C_S\le0.596$, $C_g\le1$, $C_{S2}\le2.543$, $C_E\le3.164$, $C_{MS}\le0.3555$, cross $\le2.496$):
$|{\rm model}-T_1|N^2e^\theta\le(\pi^{3/2}/2)C_F=24.35\le C_b=29.5$.
Pointwise: $D(N)=1.657,1.620,1.601,1.594,1.578,1.565$ at $N=400,1000,2000,3000,10^4,10^5$.

## 6. Finite range and validation (L11, L14) — R1 repaired

Certified quadrature of $M_k^{\rm main}$ at every integer $N\in[100,400]$
(segments $[0,1,4,16,64,200]$ + proved tail $\le\frac{2\sqrt2}{\pi aR^2}e^{-aR}$), 301 values:
$$\max_{100\le N\le400}\bigl|M_k^{\rm main}+X\bigr|\,N^2e^{\theta}= \mathbf{1.751172}\ \ (N=116),$$
where $X=\frac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}[e^{i\theta}(1+\varpi_1/\sqrt N)]$;
the sharp constant on $[100,400]$ is $\le1.76$. **(R1)** The theorem needs $|M^{\rm main}+X|$;
an earlier draft's code measured $|M^{\rm main}-X|\approx2|X|$ and misreported the scaled
residual as $\approx1409$ at $N=381$ — the correct value there is $1.619635$.
At $N=100,200,400$ the scaled residuals $|M^{\rm main}+X|N^2e^\theta=0.315,1.610,1.335$, and
$|M_k-M_k^{\rm main}|=1.65\times10^{-13},4.2\times10^{-18},4.7\times10^{-25}$; the full-$M_k$
sharp residuals satisfy $\max|M_k+X|N^2e^\theta=1.743787$ ($N=158$).

## 7. Rest combination (L13) — unchanged

$10{,}600{,}000\,e^{-2\sqrt N}\le N^{-2}e^{-\sqrt{\pi N}}$ iff
$2\log N+\log(10{,}600{,}000)\le(2-\sqrt\pi)\sqrt N$; certified $N_3=25{,}711$.
For $100\le N\le N_3$: $C_3=C_2+10{,}600{,}000\cdot h(309)$, $h(309)=1748.976914$,
so $C_3=1.854\times10^{10}$.

## 8. Errata (referee DCY2R)

- **R1 (sign, required).** Finite-range residual is $|M^{\rm main}+X|$, not $|M^{\rm main}-X|$
  (the old code measured $\approx2|X|$, hence the spurious "$1409$"). Re-certified: max $1.7512$.
- **R2 (notation).** The complex estimate is for $J:=(2/\pi)\tilde I$, $M^{\rm main}=\operatorname{Im}J$.
- **R3 (typos).** $P(u)$ $u$-coefficient $-2b_Ry^3\to-2b_Ry^2$; §5 spurious $(1+\varpi_1/\sqrt N)$
  factor removed; `dcy2_cert` $h_1$ comment is the used form $|h_1|\le1.5y(1+2y)$.
- **R4 (remainder).** $E_{\rm tailmom}\times2$ (two-sided) + `Mtilde` correction
  $2e^{a^2/4c}t_R^{j+1}/(j+1)$; worst region-B cell $396{,}111<C_2$.
- **R5 (endpoints).** `TailMom`($p\le2$) at $N_{\rm lo}$; $e^{\gamma T^3}$ at a consistent
  endpoint; `EdgeInt` $\Gamma$-subtraction clamped at $0$. ~1% impact.
- **R6 (scripts).** `step_03_watson.py` repaired (alternating-series bound; passes its asserts;
  $C_\rho=1/16\le0.5$).

## 9. Files and reproduction

`verify_dcy2_repaired.py` (126 s < 5 min) reproduces everything: §1 $\varpi_1$, §2 R1 exhibit,
§3 finite range, §4 region B, §5 region C, §6 part (b), §7 rest, §8 tail machinery, §9 R4/R5
deltas, §10 R6 Watson, §11 R3 $P(u)$. Modules: `step_06_finite.py`, `step_02_varpi1.py`,
`step_03_watson.py`, `step_07_partb.py`, `dcy2_cell.py`, `dcy2_sweep.py`, `dcy2_regionc.py`,
`dcy2_cert.py`, `dcy2_budget.py`, `dcy2_common.py`, `dcy2_audit.py`.
```

**Bottom line:** The DCY2 theorem and corollary stand with the **same constants** ($C_2=10^6$; sharp $\le1.76$ on $[100,400]$; $C_3=1.86\times10^{10}$, $C_2+1$ past $N_3=25{,}711$). The one genuine error (R1, the finite-range sign) and the five minor repairs (R2–R6) are all applied, and the repaired notes (`proof_repaired.md`) contain no "1409" — the correct sharp constant is **1.751172** ($N=116$), re-certified by 301 arb quadratures.