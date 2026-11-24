
import os, datetime
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
result = r'''# SRF3: the superellipsoid's octahedron end — radius zero PROVED, exact Borel
# singularities, and the growth law

**Work dir:** `research/08_theory/SRF3_octahedron_radius_proof/`. Python:
`.venv/bin/python` (mpmath 1.3.0, sympy 1.14, python-flint 0.9.0). Every number below was
printed by my code (step files, or `verify_srf3.py` which re-runs the essentials in 117 s).
Main writeup with the full proofs: **`proof.md`** (with a lemma ledger).

**Notation.** $p=1+d$, $\tau=1-\frac1p=\frac d{1+d}$, $S(p)$ = area of $|x|^p+|y|^p+|z|^p=1$,
$S(1)=4\sqrt3$. $z_0=\frac{-\log2+i\pi}2$, $R=|z_0|=\frac{\sqrt{\pi^2+\log^22}}2
=1.6085752558559\ldots$, $\theta=\arg z_0=1.7879528081463\ldots$

---

## Results in one place

| # | statement | status |
|---|---|---|
| A | $S(p)=8(1-\tau)^2T(\tau)$, $T(\tau)=\iint_\Delta(t_1t_2t_3)^{-\tau}\sqrt{\sum t_i^{2\tau}}\,dt_1dt_2$; equivalent square / Laplace / softmax forms | **PROVED** (+40-digit check, 6 values of $d$) |
| B | Borel transform of $T$ (and of $S/S(1)$) is analytic in $|z|<R$; its only nearest singularities are $z_0=\frac{-\log2\pm i\pi}2$, **square-root** branch points. $\log2$ = log of the edge-midpoint ratio $\frac12$ (an edge of the face-simplex has exactly 2 surviving barycentric coordinates; $1+2e^{2z}=0$) | **PROVED** |
| C | $\limsup|a_k/k!|^{1/k}=1/R=2/\sqrt{\pi^2+\log^22}=0.6216681478591508\ldots$ for $a_k=[\tau^k](S/S(1))$ (same in $d$): **radius zero** | **PROVED** (upper bound §3.2, kernel non-cancellation §3.3 à R2, lower bound §3.4) |
| D | growth law: $a_k/k!\sim-\sqrt{\frac{3\pi}8}\,k^{-3/2}\,\mathrm{Re}[z_0^{1/2-k}]$; in $d$: $a_k(d)/k!\sim+\sqrt{\frac{3\pi}4}\,k^{-3/2}\,\mathrm{Im}[z_0^{1/2-k}]$; i.e. $a_k\sim C\,k!\,A^k\,k^{\beta}\cdot(\cos,\sin)$ with $A=1/R$, $\beta=-3/2$, $C_\tau=\sqrt{3\pi R/8}=1.3766110843131\ldots$, $C_d=\sqrt2\,C_\tau$ | **PROVED** (form), **NUMERICAL** (8–10 digits, vs paper's 2–4) |
| E | Borel summable along $\mathbb R_+$ ($\mathrm{Re}\,z_0<0$): Padé–Borel resummation vs direct quadrature | **NUMERICAL** |

## The growth law vs the paper's Table II row

The paper writes $a_k\sim C\,k!\,A^k\,k^\beta$ with $A=2/\sqrt{\pi^2+\log^22}=0.621668\ldots$
and $z_0=(−\log2±i\pi)/2$, status "Numerical; 2 digits (Hankel), 4 (fit)". Now:
- $A$, $z_0$: **proved** (Theorem B/C), with the exact geometric mechanism (edge strips of the
  face-simplex; $1+2e^{2z}=0$ at the edge-midpoint ray).
- $\beta=-3/2$: proved (square-root branch) and fitted ($-1.49999\pm10^{-4}$, K=120 data).
- $C$: **closed form** $\sqrt{3\pi R/8}=\frac{(3\pi)^{1/2}(\pi^2+\log^22)^{1/4}}4
  =1.3766110843131\ldots$ ($\tau$-series), fitted to **7.9 digits** (rel err $1.4\cdot10^{-8}$);
  $d$-series $C_d=\sqrt{3\pi R/4}$ fitted to **9.7 digits**; phases
  $\cos(k\theta+\pi-\theta/2)$ ($\tau$), $\sin((k-\frac12)\theta)$ ($d$) to $10^{-8}$.
  The $B(3/2,3/2)=\pi/8$ of the simplex edge (transverse beta integral) is inside $C$.

## Evidence table

| claim | evidence (residual / agreement) | digits | script |
|---|---|---|---|
| Thm A representation | 6 values of $d$, worst 40.2 digits vs H1 direct quadrature | 40 | step_01, verify_srf3 V1 |
| Jacobians ($p^{-2}(t_1t_2)^{1/p-1}$, $1-u$, $t_1t_2t_3$) | sympy exact, 0 | exact | verify_srf3 V1b |
| zeros $\tau^*(s,u)\sim z_0/x$ | ratios $R/(x|\tau^*|)=0.85,0.90,0.92\to1$ at $s=1/2$ | — | step_06, verify V2 |
| corner fiber $\min|z(w)|$ | $\ge R$ (edges), $=2R$ at $w=\tfrac12$ | — | step_06b, verify V2b |
| kernel bound $|k_j|\le(Ce^{-V})^j/j!$ | $C(V)=2e^V|\log(1-e^{-V})|\le4.2\to2$ | exact | step_06b, verify V3 |
| engine K=120 vs paper thm:oct $d^2..d^4$ | abs err $10^{-60}$; $w_2..w_6$ exact match | 60 | step files, verify V4 |
| $\beta=-3/2$ | grid optimum $-1.49999$, windows 40–90; rms($-1.5$) $1.2\cdot10^{-10}$ vs $2.7\cdot10^{-8}$ ($-1.49$) | 4–5 | step_04, fits |
| $C_\tau=\sqrt{3\pi R/8}$ | fitted 1.3766111033561893, rel err $1.4\cdot10^{-8}$ | **7.9** | fits (K=120, nc=5, k0=90) |
| $\varphi_\tau=\pi-\theta/2$ | fitted $-2.24761626377706$, abs err $1.4\cdot10^{-8}$ | **8** | fits |
| $C_d=\sqrt{3\pi R/4}$ | fitted 1.9468220658997473, rel err $1.8\cdot10^{-10}$ | **9.7** | fits |
| $\varphi_d=\varphi_\tau+3\pi/2$ (factor $e^{-\bar z_0}=\sqrt2i$) | fitted 2.46477271476253 | 8 | fits |
| edge-only content | vertex angle 1.3536 amplitude $<2\cdot10^{-3}$ of edge | — | step_02/04 |
| model vs coefficients | max rel resid $1.9\cdot10^{-11}$ ($\tau$), $3.2\cdot10^{-12}$ ($d$), k=90..120 | — | fits |
| mixing $\tau\to d$ | exact identity $c_k=\sum_j a_j(-1)^{k-j}\binom{k-1}{j-1}$, resid $10^{-51}$ | 51 | step (mixing) |
| Padé–Borel summability | digits at $d$: 41.4 (0.001), 40.6 (0.003), 37.0 (0.01), 15.2 (0.03), 6.3 (0.1), 37.1 ($-0.01$), 13.7 ($-0.03$), 3.9 ($-0.1$) | 6–41 | step_05, verify V6 |
| 2D control | $C=\sqrt2$, $\beta=-3/2$, radius $2/\pi$, phase $-3\pi/4$ reproduced | 3–4 | step (control) |

## Coefficient artifacts
- `SRF3_hi_coeffs_K120.json`: $[\tau^k](S/S(1))$ and $[d^k](S/S(1))$, $k\le120$, ~460 dps
  (H1 Route-A engine, local copy `engine_h1.py`; cross-checked vs stored K=60 exactly).
- `w_series_K120.json`: pure-weight log-series coefficients $w_k$ ($w_2..w_6$ = paper's
  thm:oct values, verified: $-\pi^2/9,-16\zeta(3)/9,-37\pi^4/810,-160\zeta(5)/27,
  16\zeta(3)^2/81-572\pi^6/76545$).

## Step files
- `step_01.py` exact representation + 40-digit check.
- `step_02.py` frequency scan (edge vs vertex). `step_03.py` Padé–Borel poles (noisy).
- `step_04.py` annihilator + two-frequency LSQ: $\beta=-3/2$, edge-only.
- `step_05.py` Padé–Borel resummation (Task 4).
- `step_06.py` zeros $\tau^*(s,u)$. `step_06b.py` corner fiber min + K-factor bound.
- `engine_h1.py` K=120 coefficient engine (H1 Route A, local). Fit steps (inline, transcript):
  C/φ/β fits, mixing checks, law tables.
- `verify_srf3.py` consolidated verification, 117 s, all PASS.

## Relationship to prior work
- Uses (read-only): paper `main.tex` (sec:octahedron, thm:oct, thm:zero), SRF1 (numerical
  radius-zero, K=60), H1 (K=80 coefficients, direct quadrature), D16 (pure-weight engine),
  R2 (the proved 2D rhombus radius-zero, whose $K-1$ bound §3.3 mirrors), GRW1 (parallel 2D
  Borel-summability; not duplicated — only the 2D control fit was run here).
- Upgrades the paper's Table II octahedron row from "Numerical, 2–4 digits" to **PROVED with
  closed forms and 8–10 digits**.
'''
open(wd+"/result.md","w").write(result)
print("result.md written,", len(result), "chars")

# update progress.log
log = open(wd+"/progress.log","a")
more = [
 "step_05: Pade-Borel summability along R+ vs direct quadrature: 41.4 digits (d=0.001) .. 6.3 (0.1), 37.1 (-0.01) [NUMERICAL, Task 4]",
 "step_06/06b: zeros tau*~z0/x(1+O(1/x)); corner fiber min |z|>=R (edges), 2R (midpoint); K-factor bound |k_j|<=(C(V)e^-V)^j/j!, C(V)<=4.2 [NUMERICAL/exact]",
 "growth law FINAL: beta=-3/2, C_tau=sqrt(3 pi R/8) to 7.9 digits, C_d=sqrt(3 pi R/4) to 9.7, phases to 1e-8; beta=-1.49999+-1e-4 [NUMERICAL]",
 "w-series (log): identical leading law as S/S(1) (C_w/C_tau=1.00000001), w_2..w_6 = paper thm:oct exact [NUMERICAL]",
 "proof.md + result.md written; verify_srf3.py all PASS in 117s",
]
for e in more:
    log.write(f"[{datetime.date.today()}] {e}\n")
log.close()
print(open(wd+"/progress.log").read()[-800:])
