# HEX5 — Verdict: **CONFIRMED**

For the regular hexagon $K_p=\{\sum_{i=1}^3|\langle n_i,x\rangle|^p\le1\}$, $n_i=(\cos(i\pi/3),\sin(i\pi/3))$, the independent numeric computation gives
$$\boxed{\ [\zeta(5)]\,P^{(5)}(1)=\frac{2835\sqrt3}{32},\qquad [\zeta(2)\zeta(3)]\,P^{(5)}(1)=\frac{1395\sqrt3}{32}\ }$$
exactly as claimed. Both coefficients were recovered *purely numerically* (no input from the paper's edge law): a stable 18-dimensional integer relation for $P^{(5)}(1)/\sqrt3$ on the prescribed weight ≤ 5 basis, found by PSLQ with the **identical integer vector at two independent precisions (125 and 145 digits)**.

## 1. Sanity: P′′′(1) reproduced
Engine (Route A of `../CX3k_review_draft8/polygon_engine.py`, dps=60) vs the paper's closed form $(\sqrt3/8)(16\log^32+3\pi^2\log2-96\log^22+96\log2-9\zeta(3)-6\pi^2)$:

| | value |
|---|---|
| engine | −5.145455017135615837171681131123975071115553566419741**699** |
| paper | −5.145455017135615837171681131123975071115553566419741**721** |

abs diff **2.24e-53** (engine last-level diff 3.35e-56). Also P(1)=3.4641016…=2√3 = hexagon perimeter ✓.

## 2. P⁽⁵⁾(1) to ≥ 80 digits (two quadrature levels)
The integrand series is π/3-periodic and wedge-midpoint symmetric (verified to 1e-40), so one wedge ×6 with paired evaluation (optimized engine, series order N=7, exact in d). Convergence is doubly exponential: level diffs 1.3e-18 → 2.1e-41 → 6.7e-88 → {5.3e-115 (dps=115 run) / 7.2e-151 (dps=150 run)}.

```
P^(5)(1) = 656.41536196762475447664760905831120964826426464898548604185355096214310290964572236224057164257468405023297944386495...
P^(5)(1)/√3 = 378.98158593221378613681191968953666478344947302239726283776139955442837890161391759674454467680213036339689469261...
```
The two independent runs share **113 leading digits**. Cross-check with the unmodified original engine at dps=50: 656.415361968… (agrees to its full ~10-digit accuracy).

## 3. PSLQ of P⁽⁵⁾(1)/√3 (basis: 1, log2, π², log²2, ζ3, π²log2, log³2, π⁴, ζ3log2, π²log²2, log⁴2, ζ5, ζ2ζ3, π⁴log2, ζ3log²2, π²log³2, log⁵2)

| precision | outcome | residual | maxcoeff |
|---|---|---|---|
| 100 | spurious, rejected (residual/maxcoeff = 7.6e-78 ≫ maxcoeff·10⁻¹⁰⁰) | 2.74e-73 | 36116 |
| 125 | **relation R** | 7.34e-122 | 61440 |
| 145 | **R again, identical vector** | 4.81e-142 | 61440 |

Residuals pass the acceptance test (~maxcoeff·10^(−digits)); independent recheck at dps=170: residual **9.59e-146**, and R reconstructs P⁽⁵⁾(1)/√3 to **149.7 relative digits**. Li₄(1/2), Li₅(1/2) not needed.

$$\begin{aligned}128\,\frac{P^{(5)}(1)}{\sqrt3}=\;&30720\log2-4920\pi^2-61440\log^22+6240\pi^2\log2+30720\log^32+448\pi^4\\&+5760\,\zeta_3\log2-1920\pi^2\log^22-5120\log^42+\mathbf{11340\,\zeta_5+5580\,\zeta_2\zeta_3}\\&-123\pi^4\log2-1440\,\zeta_3\log^22+160\pi^2\log^32+256\log^52\end{aligned}$$
(coefficients of 1 and ζ₃ are exactly 0). Since 11340/128 = 2835/32 and 5580/128 = 1395/32 (both verified `True` in code), the claimed coefficients match **exactly**. Subtraction variant (fixing the claimed ζ-terms, PSLQ on the remaining 15 elements) returns exactly the restriction of R at dps=125.

## Files that produced this (work dir `HEX5_hexagon_p5/`)
- `step_01.py` — sanity P′′′(1) original engine vs paper (dps=60)
- `step_03.py` — π/3-periodicity + midpoint-symmetry verification
- `step_04.py` — optimized engine (one wedge ×6, paired eval, per-level log); P⁽⁰…⁵⁾(1) at dps=50/115/150
- `step_05.py` — original-vs-optimized head-to-head diagnostic (exposed a 15-digit √3/2 init pitfall, fixed)
- `step_06.py` — PSLQ at dps=100/125/145
- `step_07.py` — dps=170 residual + reconstruction + subtraction-variant PSLQ
- `step_08.py` — unmodified original engine cross-check (dps=50)
- `hex5_report.py`, `hex5_finalcheck.py` — wrote `result.md` (7404 chars) and the final consolidated check
- `p5_run115.txt`, `p5_run150.txt` — raw outputs; (`step_02.py` superseded driver; `step_09…28.py` are the tool's auto-saved transcript blocks)