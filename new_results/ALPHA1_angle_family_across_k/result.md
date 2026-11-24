# ALPHA1: the corner-angle family \(M_k(\alpha)\) across \(k\)

**Workdir:** `new_results/ALPHA1_angle_family_across_k/`. All steps run under `nice -n 10` with the venv python; each line of `progress.log` records a finished step.

## Verdict

- **Task A (general-\(c\) integrand + reducer + controls): DONE — PROVED/EXACT.**
- **Task B (polynomial law; weight two; weight three): DONE** — polynomial law EXACT at 4 rational angles (\(d\le5\), \(k\le14\)); weight-two kernel universal with closed form; weight-three slope a closed form in \(c\).
- **Task C (leading law in the angle): PARTIAL** — top Newton table extracted EXACTLY at 4 angles + squircle for \(d\le6\) (squircle \(d\le8\)); the master factor identified at \(d=2\); the general \(c\)-dependence and Mellin/Gegenbauer generating function remain CONJECTURE (data in hand, not pinned).

## 1. The general-angle setup (Task A, part 1) — PROVED algebra, NUMERICAL 50 digits (`step_01.py`)

With \(t=g/(1+g)\): \(M(u;\alpha)=\int_0^1\Psi_u\,dt\), \(\Psi_u=\frac12\big(\sqrt{t^{2u-2}+(1-t)^{2u-2}-2\cos\alpha\,(t(1-t))^{u-1}}-t^{u-1}-(1-t)^{u-1}\big)\) (at \(\alpha=\pi/2\): the squircle). With \(c=\cos(\alpha/2)\), \(r=\frac{1-c}{1+c}=\tan^2(\alpha/4)\):
\[t=\tfrac{(w-r)(w+1)}{2(1-r)w},\quad \tfrac{dt}{dw}=\tfrac{w^2+r}{2(1-r)w^2},\quad \sqrt{t^2+(1-t)^2-2\cos\alpha\,t(1-t)}=\sqrt{1-4c^2t(1-t)}=\tfrac{w^2+r}{(1+r)w}\]
so \(M(u;\alpha)=\int_r^1\Psi_u(t(w))\,dt/dw\,dw\), everything rational over \(\mathbb Q(c)\). Letters \(\{0,\pm1,\pm r\}\); the only non-real poles at \(w=\pm i\sqrt r=:\pm\theta\), \(\theta^2=-r\); involution \(w\to r/w\). Checked vs `corner_series.py`: function level and \(M_k\), \(k\le4\), 2 angles, ~50 digits.

## 2. Integrand and triangular IBP system (Task A, part 2) — EXACT, controls PASS

`[eps^k]\Psi=\sum_m\tilde R_m(t)x^my^{k-m}\) with the three-block \(u=\alpha X+\beta Y+\gamma Z\), \(\gamma=-2\cos\alpha\,t(1-t)/q_{2\alpha}\) (`cwintegrand.py`; at \(\cos\alpha=0\) it reduces to the squircle's two-block \(u\)). In \(w\): \(R_m(w)=C_k\tilde T_m(w)/[w(w^2+r)^{2k-2}(w-r)(w+1)(1-w)(w+r)]\), \(C_k=1/(2(2c)^{2k-1})\) — simple real poles at the letters, pole order \(\le2k-2\) at \(\theta\).

The analogue of `eq:triangular` (`cwred.py`, mirrors WTS `fastred`): \(W_{ij}=\varepsilon^dR_{ij}-(i{+}1)G_{i+1,j}\omega_1-[i{>}j](j{+}1)G_{i,j+1}\omega_2-[j{=}i{-}1]\,i\,\mathcal F(G_{ii}\omega_1)\), \(G=\mathcal JH(W)\) normalized \(G(r)=0\), fold \(\mathcal F(R)=rR(r/w)/w^2\), \(\omega_{1,2}\) vanish at \(\pm\theta\). (Bug found by component tests `step_03–06`: \(d\,\mathrm{lmt}\) has pole at \(-r\), not \(+r\).)

- **CONTROL 1 (NUMERICAL, `step_02.py`):** at \(c=1/2\), \(\sum_{d=1}^k\sum_{i+j=d}\int_r^1K^{(k)}_{ij}\,\mathrm{lt}^i\mathrm{lmt}^j\,dw=M_k\) from `corner_series.py` to **49–52 digits**, \(k=1..4\). No simple \(\theta\)-pole survives in any kernel.
- **CONTROL 2 — SQUIRCLE (EXACT, `step_07.py`):** at \(c=1/\sqrt2\) (field \(\mathbb Q(\sqrt2,i)\)), my \(w\)-kernels for \(k=1..8\), mapped by \(v=aw\), equal the published squircle `kern_k.pkl` EXACTLY at every key (real residues identical; \(\theta\)-coefficients scale \(a^{1-j}\)). **8/8 PASS.**

## 3. Polynomial law (Task B, part 1) — EXACT at 4 angles (`step_08.py`→`kern_{c12,c35,c13,c45}.pkl`, `step_09.py`, `step_10.py`)

Kernels computed exactly for \(k=1..14\) at \(c=1/2,3/5,1/3,4/5\) (293 s). For every key \((i,j)\), \(d=i+j\le5\), every real-pole residue: the \(d\)-th finite difference vanishes for \(k\ge d+1\) (degree \(\le d-1\)), \(\ge4\) spare points — **0 failures in 320 sequences**. Exact degree \(d-1\) attained per degree; \(k=d\) generically off the polynomial.

## 4. Weight two — DONE (`step_10,20,25`; cross-check `step_11`)

For every \(k\ge2\), \(K^{(k)}_{10}(w)=\frac1{2w}+\frac1{2(w-1)}-\frac1{w+1}-\frac1{2(w+r)}\): the **universal** (angle-independent) kernel = the squircle's \(R_{w2}\) (EXACT at 4 angles, all \(k\le14\); equivalently \((1-\varepsilon)(W_{10}-\varepsilon R_{10})=\varepsilon^2R_{w2}\)). Hence \(W_2(k;\alpha)=T(\alpha)\) constant in \(k\). Value (NUMERICAL: LLL/LS fit at 16–19 angles, 76–88 digits, holdouts clean; recovers squircle \(T\)):
\[\boxed{\,T(\alpha)=\tfrac{\pi^2}{12}-\mathrm{Li}_2(c)-\mathrm{Li}_2(-c)=\tfrac{\pi^2}{12}-\tfrac12\mathrm{Li}_2(c^2)\,}\]
independently confirmed by G1's Theorem C (\(W_2=\frac12\zeta(2)-\frac12\mathrm{Li}_2(c^2)\)).

## 5. Weight three — DONE (`step_22,30`; slope value `dbg4/dbg5`)

For \(k\ge3\) the degree-2 kernels are affine, slope kernels universal \(\times\,s(\alpha)\): \(S_{20}=s\big(\frac1{w+1}-\frac1{w+r}\big)\), \(S_{11}=s\big(\frac1w-\frac2{w+1}\big)\), \(s(\alpha)=\frac{1-c^2}{2c^2}=\frac12\tan^2(\alpha/2)\) (EXACT at 4 angles, \(k=3..14\); at \(c=1/\sqrt2\), \(s=1/2\) = squircle slope kernels). So \(W_3(k;\alpha)=W_3(3;\alpha)+(k-3)\,\mathrm{slope}(\alpha)\), \(\mathrm{slope}=s\cdot J(c)\) (NUMERICAL, **61 digits** at 6 angles incl. squircle control; EXACT 6-term HPL form via a04's GPL→HPL machine, `step_30.py`):
\[\boxed{\,J(c)=\log(1-c^2)\big[\mathrm{Li}_2(c^2)-\zeta(2)\big]+2\mathrm{Li}_3(1-c^2)-2\zeta(3)+\log(c^2)\log^2(1-c^2)=2T(\alpha)\log(1-c^2)-\!\int_0^{c^2}\!\tfrac{\log^2(1-u)}udu\,}\]
At \(c=1/\sqrt2\): \(-(\zeta(3)/8+\pi^2\log2/24+\log^32/12)\), the squircle slope. So the slope **is** a closed form in \(c\) (\(\mathrm{Li}_3(1-c^2),\mathrm{Li}_2(c^2)\), logs).

## 6. Leading law (Task C) — PARTIAL (`step_31.py`→`beta_tables.pkl`, `step_32.py`→`beta_squircle.pkl`)

Top Newton table \(\beta_d^{(ij)}(\rho)=[k^{d-1}]\,\mathrm{res}_\rho K^{(k)}_{ij}\) extracted EXACTLY at the 4 rational angles (\(d\le6\)) and the squircle (\(d\le8\), from `kern_k.pkl`, \(k\le18\)). Findings: within each \((c,d)\) the table is the squircle's binomial\(\times\)zigzag pattern on poles \(\{0,\pm1,\pm r\}\) times a \(c\)-dependent master factor \(F_d(c)\): e.g. \(F_2(c)=\frac{1-c^2}{c^2}=\tan^2(\alpha/2)\) (EXACT, all 4 angles); \(F_3\): \(5/2,19/12,25/4,43/64\) at \(c=1/2,3/5,1/3,4/5\) (exact data, closed form not yet identified); \(d\ge4\) tables in hand (`beta_tables.pkl`). The value-level assembly \(c(w;\alpha)=\frac1{(d-1)!}\sum\beta\,N_{(i,j),\rho}\) and the Mellin/Gegenbauer generating function are **not yet pinned** — this is the open piece (the machinery to finish it is in place).

## Files
- Modules: `cwfield.py` (fields \(\mathbb Q(\theta)\), QSI), `cwintegrand.py` (3-block \(T_m\)), `cwred.py` (reducer), `cwbuild.py`, `lllkernel.py`, `hpl_tools.py` (copied), `corner_series.py` (copied).
- Steps: `step_01` (setup, 50 d), `step_02` (control vs corner_series), `step_03–06` (component tests/DL2 bug), `step_07` (squircle control PASS), `step_08` (kernels \(k\le14\) ×4 angles), `step_09/10` (polynomial law), `step_11–20` (weight two), `step_22–30` (weight three), `step_31/32` (\(\beta\) tables). Data: `kern_*.pkl`, `beta_*.pkl`.
- `verify_alpha1.py` as a single <15 min PASS/FAIL runner is **not yet packaged**; the controls it should assert are exactly `step_01` (50-digit integrand), `step_02` (weight sums = \(M_k\), 4 levels), `step_07` (squircle kernels, 8 levels), which together rerun in ~2 min.

## Evidence table (selected)

| Claim | Type | Evidence |
|---|---|---|
| \(w\)-integrand = corner \(M(u;\alpha)\) | NUMERICAL | 2 angles, \(M\) and \(M_k\le4\), 50 digits (`step_01`) |
| Reducer weight sums \(=M_k\) | NUMERICAL | \(c=1/2\), \(k=1..4\), 49–52 digits (`step_02`) |
| Squircle control | EXACT | \(k=1..8\), all keys, residual 0 (`step_07`) |
| Polynomial law \(\deg=d-1\) | EXACT | 4 angles, \(d\le5\), \(k\le14\), 0/320 failures (`step_09`) |
| \(W_2(k)=T(\alpha)\), kernel universal | EXACT (kernel) + NUMERICAL (value) | 4 angles, \(k\le14\); \(T=\frac{\pi^2}{12}-\frac12\mathrm{Li}_2(c^2)\) at 19 angles, 76–88 digits (`step_20`) |
| \(W_3\) affine; slope \(=s(\alpha)J(c)\) | EXACT (kernel) + NUMERICAL (value) | 4 angles, \(k\le14\); \(J\) closed form 61 digits ×6 (`step_22,30`) |
| \(\beta_{d}\) tables, \(F_2(c)=\tan^2(\alpha/2)\) | EXACT | 4 angles \(d\le6\), squircle \(d\le8\) (`step_31/32`) |
| Mellin/Gegenbauer leading law | CONJECTURE | untested |