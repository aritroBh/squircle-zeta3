# DCY2R referee verdict: **PROVED AFTER REPAIRS**

The DCY2 theorem (main term, $C_2=10^6$) and its corollary (full $M_k$, $C_3=1.854\times10^{10}$) stand. One lemma (**L11**, the finite-range check) was **WRONG as coded** — a sign error — and is repaired and re-certified here (true constant ≤ **1.76**). All other lemmas are CORRECT or need only minor repairs; nothing threatens the constants.

## The discrepancy (planner finding) — RESOLVED

DCY2's finite-range code (`step_06_finite.py`, `verify_dcy2.py` §3) computes `res = abs(M - (-T1))` with `T1 = -X`, $X=(\sqrt\pi/N)e^{-\theta}\mathrm{Re}[e^{i\theta}(1+\varpi_1/\sqrt N)]$ — i.e. $|M_k^{main}-X|$ where the theorem needs $|M_k^{main}+X|$. Since $M_k^{main}\approx -X$, DCY2's `res` $\approx 2|X|$, whose scaled value $2\sqrt\pi N|\mathrm{Re}[\cdots]|$ grows like $N$, maximized at the top of the range: $2\sqrt\pi\cdot381\cdot|\mathrm{Re}[e^{i\theta}(1+\varpi_1/\sqrt{381})]|=1409.47$ — exactly their "constant ≈ 1409" / "res/bound 0.00141". The planner's 1.744 is the true quantity with the full $M_k$. Printed (arb, outward, N=381):

| quantity | value |
|---|---|
| $\lvert M_{381}^{main}-X\rvert N^2e^{\theta}$ (what DCY2 measured) | **1409.4657** |
| $\lvert M_{381}^{main}+X\rvert N^2e^{\theta}$ (the theorem) | **1.6196351** |

As committed, L11 certifies only the trivial $2|X|\le C_2N^{-2}e^{-\theta}$; $N\in[100,400]$ was **not** covered. L14's quoted residuals (0.000165, 0.000751, 0.000932) are the same artifact.

## Lemma verdicts

| # | Lemma | Verdict | Key evidence (all printed by code) |
|---|---|---|---|
| 1 | Ray representation (CX6) | CORRECT | code integrand matches $H E Q_k\,ds/s$ exactly; table cross-check $|M_k-M^{main}|\sim10^{-13}$ |
| 2 | $\varpi_1$ closed form | CORRECT | independently re-derived ($G_1,G_2,c_1,S_1$ from $\phi$; diff vs paper = 0 exactly; 40-digit match) |
| 3 | Factorization, modulus, single-max, tail rates | MINOR FIX | factorization exact but for $J:=(2/\pi)\tilde I$ ($M^{main}=\mathrm{Im}\,J$); §1 writes $\tilde I$ for both ($\tilde I\approx(\pi/2)T_1$, so " $|\tilde I-T_1|\le C_2N^{-2}e^{-\theta}$ " as printed is false). $P(u)$ $u$-coefficient is $\sqrt2(N-2)y^2-2b_R{\bf y^2}$ not $-2b_Ry^3$ (same typo in `lamL` code; harmless — min selects $\pi y/2$). $P(1)=-y(4+8\sqrt2y+\cdots)/(2Q_0)<0$ exact; $P''<0$; rates rDelta (0.514), eLamD (1.117), $\lambda_L$ verified (grids + arb) |
| 4 | $\ell$-boxes, $c_2$, $\mathrm{Re}\,\ell_1$ | CORRECT | $\mathrm{Re}\,\ell_1=-y(1+b_Ry/Q_0)$ exact (symbolic); $c_2\ge0.9\,b_Ry$ (min ratio 0.9647); $|\ell_m|/L_m\le0.994$ for $m=2..7$ |
| 5 | Watson constants | GAP (script) / constants OK | `step_03_watson.py` **fails its own assert**: Cauchy gives $M/4=0.9151>0.5$; interval part returns 0 (misuses `arb(lo,hi)` as (mid,rad)). Truth: $\max|\rho_3|/u^3=0.0117$ on $(0,1]$, $0.053$ on $[1,20]$ — so $C_\rho\le0.5$, $\le2.42$ true and easily re-certified. Used constants certified here: $C_G$: truth 0.746/0.768/0.779 ≤ 1; $|h_j|/(\tfrac32)_jy^j=0.978/0.957/0.936\le1+2y$; $M_4'$ spot-checks ≤ 0.045 |
| 6 | Moments; model $=T_1+O(\cdot)$; the two exact identities | CORRECT | both identities re-verified exactly (diff=0); $D(N)$ pointwise 1.657→1.565; moment-sum $O(y)$ coefficient = 0 |
| 7 | Remainder $E(t)$ | GAP (repairable) | kept/dropped monomial algebra verified symbolically ($P$, drop, $\beta_2,\beta_3$ consistent); audit $E_{mom}\ge$ truth ratios 0.001–0.006. But `TailMom` uses the **one-sided** formula $\frac12(c_2^m/2)^{-s}\Gamma(s,z)$ for a two-sided integral (valid bound $=2\times$); `Mtilde` misses the provable correction $2e^{a^2/4c}t_R^{j+1}/(j+1)$ ($j=0$ margin 3.5%, unused; $j\ge4$ margin ≥30%). Repair: $E_{tailmom}\times2$ + correction ⇒ worst pointwise LHS 320,253→325,886 (+1.8%); worst cell 390,055→≈397k. $C_2$ unchanged |
| 8 | Window-tail + single-dip ($W$) | CORRECT | $W(1)<0$, $W''>0$ logic sound; rate inequalities verified directly on $[\Delta,30\rho]$ + asymptotics |
| 9 | Region B sweep | CORRECT + MINOR FIX | 279 cells reproduced: worst **390,055.4767** at (3075.23, 3136.74); coverage $[400,10^5]$ exact. ~1% endpoint gaps (harmless): `TailMom`($p\le2$) decreasing but taken at Nhi (ratios 0.983/0.990/0.997 on a ≤0.4% share of LHS); $e^{\gamma T^3}$ at Nlo while $\gamma T^3$ grows ~1.2%/cell (compensated ≥1.01 by `Mtilde` at Nhi); assert omits $+C_b=29.5$ (absorbed: 390,084.5 < $C_2$) |
| 10 | Region C | CORRECT | cell-form LHS($10^5$) = **74,115.35**; monotone on a 923-point dense grid $[10^5,10^7]$; values 26,348 / 17,592 / 14,711 / 14,090 / 14,014 at $10^7,10^9,10^{12},10^{16},10^{20}$ (pointwise form); piecewise $y$-exponents $\ge y^{3/2}$, log-derivative conditions at $\kappa^2=28.2$ verified |
| 11 | Finite range $[100,400]$ | **WRONG (as coded)** → REPAIRED | sign error (above). Re-certified (301 arb quadratures): max $|M^{main}+X|N^2e^\theta$ = **1.751172383** (N=116) |
| 12 | Part (b), $C_b=29.5$ | CORRECT | all six constants verified: $C_S$ 0.596 ✓, $C_g$ 1 (truth 0.78) ✓, $C_{S2}$ 2.543 (truth 1.135) ✓, $C_E$ 3.164 (truth 0.67) ✓, $C_{MS}$ 0.3555 (truth 0.29) ✓, cross 2.496 (truth 2.4951) ✓; $C_F=8.7457\le10.155$; $C_b=25.375\le29.5$. (§5's factor $(1+\varpi_1/\sqrt N)$ is spurious; correct: $|{\rm model}-T_1|N^2e^\theta\le(\pi^{3/2}/2)C_F$) |
| 13 | Rest combination (CX6) | CORRECT | $N_3=25{,}711$ (boundary: 1.00012 > 1 vs 0.99948 ≤ 1); $h(309)=1748.976914$; $C_3=1.8540155\times10^{10}$ |
| 14 | Validation at 100,200,400 | WRONG numbers (artifact of #11); conclusion right | true scaled residuals: 0.315, 1.610, 1.335 |

**Audit-script caveat** (`dcy2_audit.py`): its mpmath ray quadrature is unreliable for $N\gtrsim5\cdot10^3$; its "FAIL at $N=10^5$, ratio 2.046" is a numerical artifact. Reliable window-decomposition truth (cross-checked to 8+ digits at N=400/2000): $|J-{\rm model}|N^2e^\theta = 0.00092,\ 0.00055,\ 0.00027,\ 8.7{\times}10^{-5},\ 2.8{\times}10^{-5}$ at $N=400,2000,10^4,10^5,10^6$ — decreasing; certified bounds dominate truth by $10^8$–$10^9\times$.

## Repairs (all written out)

- **R1 (required, done):** fix the sign in `step_06_finite.py`/`verify_dcy2.py` §3 (`res = abs(M - T1)`). Re-certified: max 1.7512 ≤ $C_2$; full-$M_k$ table residuals max **1.743787030** (N=158), matching the planner.
- **R2 (write-up):** state the complex estimate for $J:=(2/\pi)\tilde I$, $M^{main}=\mathrm{Im}\,J$ (the code's "Ĩ/model/target" are all $J$-quantities).
- **R3 (typos):** $P(u)$: $-2b_Ry^3\to-2b_Ry^2$ (text + `lamL_of`/`Pofu`); §5 spurious $(1+\varpi_1/\sqrt N)$ factor; `dcy2_cert` comment $|h_1+\tfrac32y|\le y^2$ (used form: $|h_1|\le1.5y(1+2y)$).
- **R4 (moments):** $E_{tailmom}\times2$ (provably-valid two-sided); `Mtilde` correction $2e^{a^2/4c}t_R^{j+1}/(j+1)$ (or prove the empirical domination). Impact ≤ +2% (worst cell), +17.6% (N=400).
- **R5 (sweep endpoints):** `TailMom`($p\le2$) at Nlo; $e^{\gamma T^3}$ at consistent endpoint (or ×1.015); `EdgeInt` $\Gamma$-subtraction convention. ~1% impact.
- **R6 (scripts):** repair `step_03_watson.py` (Cauchy $M/6$ + proper interval cells); replace the audit's ray quadrature for $N\gtrsim5\cdot10^3$ by the window decomposition.

## Proposed paper statement

> **Theorem (first correction, computer-assisted).** Put $N=k+1\ge100$, $\theta=\sqrt{\pi N}$, and $\varpi_1=\bigl[(24+12\pi+\pi^2)+i(24-12\pi+\pi^2)\bigr]/(48\sqrt\pi)=0.84121508776758439708\ldots-0.04501183768517361657\ldots i$. Then, with $M_k^{\rm main}$ the $\pm i\pi/2$ branch contribution,
> $$\Bigl|M_k^{\rm main}+\tfrac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\tfrac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|\le 10^{6}\,N^{-2}e^{-\theta}\qquad(N\ge100),$$
> and, adding the proved rest bound $|M_k-M_k^{\rm main}|\le10{,}600{,}000\,e^{-2\sqrt N}$,
> $$\Bigl|M_k+\tfrac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\tfrac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|\le C_3\,N^{-2}e^{-\theta},\quad C_3=1.86\times10^{10}\ (100\le N\le25{,}711),\quad C_3=10^6+1\ (N>25{,}711).$$
> On $100\le N\le400$ the certified constant is ≤ 1.76. Proof: 301 certified quadratures; 279 interval cells on $[400,10^5]$ (worst $3.91\times10^5$); a monotone analytic tail on $[10^5,\infty)$. In particular the $O(N^{-1})$ term and the value of $\varpi_1$ in obs:Mk are now proved; only the rest exponent $\sqrt{2\pi}$ remains numerical.

## Files (cwd: `new_results/DCY2R_referee_varpi1`)

- **`verify_dcy2r.py`** — runs in **132 s** (< 5 min), all checks pass: §1 $\varpi_1$ (independent sympy derivation, exact diff = 0); §2 discrepancy exhibit (N=381, both signs); §3 repaired finite range (301 arb quadratures + certified table); §4 region B (their `LHS_cell`, 279 cells, worst 390,055.4767); §5 region C (LHS(1e5)=74,115.35, monotone grid); §6 part (b) constants ($C_F=8.7457$, $C_b=25.375$); §7 rest ($N_3$, $h(309)$, $C_3$); §8 tail machinery (arb, outward); §9 repair deltas.
- **`result.md`** — the full referee report (lemma table, repairs, certified constants).
- **`progress.log`** — one line per step.
- **`step_01.py`–`step_67.py`** — the working computations: independent $\varpi_1$ derivation (step_25); first finite-range repair (step_21); region B rerun (step_34); extended audit incl. the $N=10^5$ artifact (step_35, step_39–42); window-decomposition truth (step_44–46, step_52); certified arb $|J-{\rm model}|$ attempts (step_50); tail-rate grids (step_38, step_48); Watson/$c_2$ certification (step_51); part-(b) constants (step_55, step_60); repair-delta quantification (step_49); TailMom/endpoint issues (step_57, step_58); region C reconciliation (step_61); rest combination (step_54); M₄′ checks (step_63, step_65); kept-polynomial algebra (step_62); $P(u)$ corrected coefficient + $P(1)<0$ exact (step_66); L4 checks (step_67); verifier timing runs (step_53, step_56, step_59).

Every number quoted above was printed by code in this session.