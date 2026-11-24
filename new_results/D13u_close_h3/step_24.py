md = r"""# D13u — closing the last write-up gap of the decay theorem (H3 lower box)

**Scope.** Second D13t review finding (result.md §R5 of `D13t_review_fixes/`): the *tight* h3 lower box
of `fixed/step_38.py` contains the truth `Re Φ(x0)+θ` on the saddle curve (σ* ∈ [0.4286, 0.5)) but **fails
at σ = 1 for N ∈ [100, 130.68]**. Question: does the certificate need the box over the whole σ-subcell
(Task 1)? Does the *coded* box contain the truth everywhere anyway (Task 2)? Amend H3 (Task 3).

**Verdict.** Only the **saddle value** must lie in the mS box (Task 1). Independently, the coded box is
certified to contain the truth over the **entire** (y,σ) box of all 29 cells × 4 σ-subcells (Task 2) —
no widening needed, `step_18` sup **unchanged at 194.1267251** (verbatim rerun, max |ΔC| vs stored
`tight_100_400.json` = 0.0). The theorem constant C = 254.2276 stands. Amended H3 in §3.

All numbers below were printed by my code (`step_NN.py` in this directory; python-flint arb/acb @ 60–70 dps,
mpmath @ 50–80 dps). `fixed/` = `.../D13s_fix_recertify/D13s_work/fixed/`.

---

## Task 1 — which quantity the mS / h3 box must contain: **the value at the saddle only**

`fixed/step_38.py` builds the box (whole (y,σ) subcell) and then uses it in exactly two places:

```python
25: def h3encl(T):
26:     T = arb(T).upper()
27:     return iv(rt2/3/(1+rt2*T+T*T), rt2/3)
...
50:     rho = 1 + y*y*sig
51:     y0 = rho*y
52:     x0 = 1/y0
54:     N = pi/(2*y*y); th = pi/(rt2*y)
59:     mS = (N*y0**3/2)*h3encl(y0)
60:     eRe0 = arb(mS.exp().upper())          # = exp(RePhi0+th) upper
70:     P0t = g0m_up*eRe0
134:     TR = GR*eRe0*arb((-c2R*dlo**2).exp().upper())/arb(lamR.lower())
140:     TLm = GLm2*eRe0*arb((-c2L*dlo**2).exp().upper())/arb(lamL.lower())
182:     z = (2*N*I1/pi**arb("1.5"))*g0box*mS.exp()*((-2*pi*invp2).sqrt())*((om**2*invp2/2).exp())*eiphi
183:     rhoz = arb(abs(z-1).upper())
184:     D = sqpi*rhoz*arb(N.sqrt().upper())
186:     C = (2/pi)*E_I*N15 + D
197:     for k in range(nsub):
198:         sig = iv(arb(k)/nsub, arb(k+1)/nsub)
199:         cell(y, sig, J, c_w, coll)
```

**Use (i) — `eRe0`, line 60: only the UPPER endpoint.** `eRe0 = exp(mS).upper()` bounds
`exp(Re Φ(x0)+θ)` at the reference saddle point `x0 = 1/(ρy)` of the *true* parameters, and feeds only the
error integrals `W1, W2, W3` (via `P0t`, line 70), `TM` (line 128, via `P0t`), `TR` (134), `TLm` (140).
In each, the integrand is factored `e^{Re Φ(x)+θ} = e^{Re Φ(x0)+θ} · e^{Re Φ(x)−Re Φ(x0)}` and the second
factor is bounded by the curvature machinery (`c2, c2R, c2L, lamR, lamL`, lines 88–143) — the x-variation
never touches `mS`. So `eRe0` is needed only at the true point `(y_true, σ*)`. The lower endpoint of `mS`
is irrelevant for use (i).

**Use (ii) — `z`, line 182: the full interval, but only at the saddle.** `z` is the saddle-model ratio of
the main term to `(√π/N)e^{−θ+iθ}`; `rhoz = |z−1|.upper` (183) must dominate `|z_true−1|`, where `z_true`
is evaluated at the true saddle of the cell. Every factor of `z` is a box that must contain the *saddle*
value: `g0box ∋ g(x0)` (midpoint + G1-propagation, lines 61–66), `eiphi ∋ e^{iφ(saddle)}` (181), and
`mS.exp() ∋ e^{Re Φ(x0)+θ}` — i.e. **`mS.lower ≤ truth(y_true,σ*) ≤ mS.upper`**, both endpoints.

**Why whole-σ-box containment is not required.** The saddle satisfies σ* ∈ [0.428579581, 0.5) for all
N ∈ [100, 10²⁰] (D13t R1; reproduced: [0.428579581005, 0.460474347638] on N ∈ [100, 400], step_14), so the
true point always lies in σ-subcell k = 1 ([0.25, 0.5]) of the `nsub = 4` split (lines 197–198), and
`certified_y38` takes the max of the per-subcell `C` (line 187) — only subcell k = 1 is load-bearing; the
other subcells' `C` merely loosen the max. The σ = 1 lower-box failure (a) lives in subcell k = 3, which
carries no truth, and (b) would be harmless even in k = 1 as long as the saddle value is contained.

**Corollary.** The certificate is sound iff `truth(saddle) ∈ mS` (both endpoints) and
`truth(saddle) ≤ mS.upper` for `eRe0` — i.e. containment **at the saddle only**. R5's domain failure is
therefore immaterial; and the upper box is valid globally anyway (worst upper margin +2.0207e-3 over the
whole cover, Task 2), so `eRe0` is safe with room to spare.

---

## Task 2 — the coded mS box contains the truth over the WHOLE (y,σ) box anyway (all 29 cells × 4 subcells)

Since Task 1 says whole-σ containment is not required, no repair is needed for soundness; but R5's option
(b) — certify the coded box directly — is achievable and closes the write-up gap completely. Checks below
use the exact identity `truth(y,σ) := Re Φ(x0)+θ = −(N/2)S(t) − (π√2/(4t))(ρ−1)² = (N/2)ψ(t) − (π√2/(4t))(ρ−1)²`,
`t = ρy`, `ψ(t) = √2 t − log(1+√2 t + t²)`, with `ψ` evaluated cancellation-free as `−log(1−w)`,
`w = 1 − (1+√2 t + t²)e^{−√2 t}`.

1. **Identity check** (step_11): truth via identity vs direct `acb` evaluation of
   `Re[(−1+i) a x0 + N log(b x0/(1+b x0))] + θ` at 600 random (N,σ), N ∈ [100,400]:
   worst relative residual **3.962977836279654e-66** @ 70 dps. (Also pure algebra: `a ρ²/t − a(ρ−1)²/t = a(2ρ−1)/t`, `√2 t N/2 = aρ²/t`.)
2. **Point checks** (step_11): every cell × every σ-subcell; σ on the grid {j/16, j = 0..16}
   (⊇ {0, 0.25, 0.5, 0.75, 1}) plus **200 random σ per subcell** (800 per cell, 23 200 total), each at
   y ∈ {y1, y2} ∪ 3 random y — **118 900 arb evals @ 70 dps**:
   - worst lower margin `truth − mS.lower` = **+7.896665667541336e-4** at cell [392.013, 400], subcell k=3, σ = 1, y = y1;
   - worst upper margin `mS.upper − truth` = **+2.021149884313004e-3** at cell [392.013, 400], k=0, σ = 0.25, y = y2.
   All positive; **no failures ⇒ no widening of `h3encl` is needed.**
3. **Certified continuum coverage** (step_13): rigorous mean-value (centered) enclosures of `truth` on a
   64×64 (σ,y) mesh per subcell — **475 136 box evals, 0 boxes not contained**; worst lower margin
   **+7.891687500990234e-4** (cell [392.013, 400], k=3, σ→1, y→y1), worst upper **+2.020668607589575e-3**
   (same cell, k=0). Worst subcell refined to 256×256 (step_15): certified lower margin
   **+7.896345974019293e-4**, max truth-enclosure width **1.1523e-6** (vs corner point value +7.8966657e-4 —
   overestimation ≤ 3.1e-8). Derivative bounds validated vs central differences: worst rel residual
   **2.100908942946879e-40** (300 pts, FD floor h²=1e-40).
   (Note: step_12's naïve direct interval evaluation of `truth` produces ~1e-2-wide enclosures from
   dependency in `(1+√2 t+t²)e^{−√2 t}` and spurious negative "margins"; superseded by the mean-value form.
   Retained for transparency.)
4. **Saddle containment** (step_14; the quantity Task 1 shows is load-bearing): for every cell,
   σ*(y) solved (saddle branch ρ(t) = √(1+t²/(1+√2 t))) at y ∈ {y1, mid, y2}; σ* ∈
   **[0.428579581005, 0.460474347638] ⊂ subcell k=1** for all 29 cells; truth(σ*) ∈ mS box of subcell k=1
   with worst margins **+9.868685738423095e-4 (lower) / +2.080900996286550e-3 (upper)**. Cell [100,105]:
   margins +0.005222/+0.008836 — reproduces D13t step_39's +0.00522/+0.00884 exactly.
5. **Tight-vs-coded at σ = 1** (step_20; N = n1 corner): the tight lower box `blo = (N t³/2)(√2/3)/(1+√2 t+t²)`
   fails (margin < 0) for N < 130.68, while the coded `mS.lower` — whose only slack at the (y2, σ=1) corner is
   the factor `(y0.lower/t)³ ≈ 0.9167…0.9671` (the `N.lower` and `T = y0.upper` factors are 1 there) — contains:

| cell N | truth(y2,1) | tight_lo | tight margin | coded mS.lower | coded margin |
|---|---|---|---|---|---|
| [100, 105] | 0.040374 | 0.040647 | **−0.000274** | 0.036582 | **+0.003792** |
| [105, 110.25] | 0.039544 | 0.039755 | −0.000211 | 0.035820 | +0.003724 |
| [120, 126] | 0.037340 | 0.037410 | −0.000070 | 0.033806 | +0.003534 |
| [130, 136.5] | 0.036066 | 0.036070 | −0.000004 | 0.032646 | +0.003419 |
| [373.35, 392.01] | 0.022399 | 0.022182 | +0.000217 | 0.020341 | +0.002058 |
| [392.01, 400] | 0.021897 | 0.021683 | +0.000214 | 0.020903 | +0.000994 |

   (margins +0.003792/+0.003724/+0.003534/+0.003419 reproduce D13t R5's four cells +0.00379/+0.00372/+0.00353/+0.00342).
   The worst coded margin over the cover sits in the *last* cell — its ratio is truncated (400/392.013 = 1.0204),
   so its y-box (hence `(y0.lower/t)³` slack) is smallest.
6. **Rerun of `fixed/step_18.py`** (verbatim subprocess, step_19): sup C on [100,400] = **194.1267251** at
   (373.3456322341576, 392.0129138458655); max |ΔC| vs stored `tight_100_400.json` over 29 cells = **0.000e+00**.
   Unchanged, as expected since no widening was needed.

---

## Task 3 — amended H3 statement (drop-in replacement)

**(H3, amended).** Let `ρ = 1+y²σ`, `t = ρy`, `N = π/(2y²)`, `θ = π/(√2 y)`, `x0 = 1/t`, and
`ψ(t) := √2 t − log(1+√2 t + t²)`. Then the exact identity

  `Re Φ(x0) + θ  =  (N/2)·ψ(t)  −  (π√2/(4t))·(ρ−1)²`                                 (1)

holds for all t > 0, ρ ≥ 1 (residual ≤ 4e-66 @ 70 dps on 600 random points; also exact algebra). The coded
enclosure is `mS = (N·y0³/2)·[ (√2/3)/(1+√2 T+T²), √2/3 ]` with `T = y0.upper`, `y0 = ρy` (boxes).

* **(a) Upper box — valid on the whole domain.** `truth ≤ (N/2)ψ(t) ≤ (N/2)(√2/3)t³ ≤ mS.upper`;
  certified over all 29 cells × 4 subcells with worst margin **+2.0207e-3** (and `box_hi − truth > 0` for all
  t ∈ [10⁻¹⁰, 10¹⁰] on the saddle branch, `fixed/verify_mS.py`). Hence `eRe0 = exp(mS).upper()` is a valid
  bound for `exp(Re Φ(x0)+θ)` — the only mS fact the error integrals use.
* **(b) Lower box — proved on the saddle curve (all that is used).** On the saddle branch
  `ρ(t) = √(1+t²/(1+√2 t))` (D13t R3), the normalized margin
  `H(t) := ψ(t) − √2 t(ρ(t)−1)²/ρ(t)² − (√2/3)t³/(1+√2 t+t²)` satisfies `H(0) = 0` and
  `H′(t) = t³·Hp3(t)` (residual 2.8e-61, D13t R5), and a single arb box evaluation certifies
  **`Hp3([0, 0.13]) ≥ 0.1538363242986387 > 0`** (`fixed/certify_H3.py`). Hence `H ≥ 0` on [0, 0.13]:
  the tight h3 lower box — and a fortiori the coded `mS` box — contains `truth(saddle)` whenever the saddle
  `t ≤ 0.13`. On the certified domain N ≥ 100 the saddle `t ≤ y(100)·(1+y(100)²) = 0.1273001149747653 < 0.13`. ∎
  (Numeric saddle margin, N ∈ [100,400]: min **+2.1167e-6** tight, **+9.8687e-4** coded.)
  Since σ* ∈ [0.428579581, 0.5) ⊂ subcell k = 1 of every cell and only subcell k = 1 is load-bearing
  (Task 1), saddle containment is all the certificate needs.
* **(c) Off-saddle domain statement — false, not needed, and repaired by the coded box.** The *tight* lower
  box fails at σ = 1 exactly for N ∈ [100, 130.6848917]: `G(y,σ) := ψ(t) − (√2/3)t³/(1+√2 t+t²) − √2 t(ρ−1)²/ρ²`
  has min **−5.47229038992e-6** at (y(100), σ = 1), single sign crossing at N = 130.6848917, and
  `G(y(N),1) > 0` for all N > 130.68 (in particular the whole sweep region N ≥ 400; +4.06e-13 at N = 10⁶).
  This does not affect the certificate (saddle-only usage, Task 1); moreover the coded box's built-in slack
  (lower endpoints `N.lower·y0.lower³/2` and `T = y0.upper ≥ t`) contains the truth over the entire (y,σ)
  box of every one of the 29×4 subcells — certified by rigorous mean-value interval enclosures
  (475 136 boxes, 0 violations; worst margins **+7.8917e-4 lower / +2.0207e-3 upper**; worst subcell
  refined: **+7.8963e-4**). No change to `h3encl`; `step_18` sup remains **194.1267251**.

**Consequence for the theorem.** H3 as amended is proved; combined with D13t's R1–R4 (H1, H2, A1, B_true all
HOLD) and the exact rerun of the [100,400] cover, the certified sup values stand:
max(194.1267251 [100,400], 254.2275176 [400,10¹⁰], 47.0701250 [10¹⁰,∞)) = 254.2275176 ≤ **C = 254.2276**.

---

## Evidence table (all printed by my code)

| check | result | tool / dps | file |
|---|---|---|---|
| identity (1) vs direct Φ, 600 random pts | worst rel resid 3.962977836279654e-66 | arb/acb @ 70 | step_11 |
| point checks: 29 cells × 4 subcells × (σ-grid j/16 ∪ 200 random σ) × 5 y | 118 900 evals, **0 failures**; worst lower +7.896665667541336e-4 (cell [392.013,400], k=3, σ=1, y=y1); worst upper +2.021149884313004e-3 (k=0, σ=0.25, y=y2) | arb @ 70 | step_11 |
| mean-value mesh 64×64/subcell | 475 136 boxes, **0 not contained**; worst lower +7.891687500990234e-4; worst upper +2.020668607589575e-3 | arb @ 70 | step_13 |
| derivative bounds vs finite differences, 300 pts | worst rel resid 2.100908942946879e-40 (FD floor h²=1e-40) | arb @ 70 | step_13 |
| worst subcell refined 256×256 | certified lower margin +7.896345974019293e-4; max enclosure width 1.1523e-6 | arb @ 70 | step_15 |
| saddle containment, all 29 cells | σ* ∈ [0.428579581005, 0.460474347638] ⊂ k=1; worst margins +9.868685738423095e-4 / +2.080900996286550e-3; cell [100,105]: +0.005222/+0.008836 (= D13t step_39) | arb @ 70, mp @ 60 | step_14 |
| tight box at σ=1 vs coded box | tight fails for N<130.68 (margin −2.74e-4 at N=100); coded contains (margins +9.94e-4 … +3.79e-3); slack factor (y0.lower/t)³ = 0.9167…0.9671 | arb @ 70 | step_20 |
| G(y,σ) domain scan | min G(y,1) = −5.47229038992e-6 at N=100; G(y,1)=0 at N=130.6848917 (single crossing, 6000-pt scan); min G over σ∈[0,0.5] = +2.56497258365e-14; G(y(N),1)>0 for N=400…10⁶ | mp @ 80 | step_16, step_22 |
| tight margin on saddle curve, N∈[100,400], 2001 pts | min +2.11665171625e-6 at N=400 | mp @ 80 | step_16 |
| `fixed/certify_H3.py` rerun | Hp3([0,0.13]) lower = 0.153836324298638698272433498578394058535265579455371823206015 > 0; t_max = 0.1273001149747653 < 0.13 | arb @ 60 (verbatim) | step_16 |
| `fixed/step_18.py` verbatim rerun | sup C = 194.1267251 at (373.3456322341576, 392.0129138458655); max |ΔC| vs stored tight_100_400.json = 0.000e+00 | arb @ 70 (verbatim) | step_18, step_19 |
| naïve interval `truth` box (artifact) | ~1e-2-wide enclosures → spurious negative margins; dependency in (1+√2t+t²)e^{−√2t}; superseded by mean-value form | arb @ 70 | step_12 (transparency) |

## Rerun commands

```bash
PY=<project>/.venv/bin/python
FIX=<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed
WRK=<project>/research/08_theory/D13u_close_h3
cd $FIX && $PY certify_H3.py          # Hp3([0,0.13]) lower 0.153836... ; t_max 0.1273001149747653
cd /tmp && $PY $FIX/step_18.py        # 29 cells, sup C = 194.1267251 (writes tight_100_400.json in /tmp)
cd $WRK && $PY step_11.py             # identity residual + 118 900 point checks (0 failures)
cd $WRK && $PY step_13.py             # 475 136-box mean-value certification (0 not contained)
cd $WRK && $PY step_14.py             # saddle containment, all 29 cells
cd $WRK && $PY step_15.py             # 256x256 refinement of worst subcell
cd $WRK && $PY step_16.py             # R5 tight-box numbers + certify_H3 reproduction
cd $WRK && $PY step_20.py             # tight-vs-coded table at sigma=1 + saddle detail
cd $WRK && $PY step_22.py             # saddle margins (full digits) + G(y,1) sign structure
```

## Files produced (D13u_close_h3/)

`step_11.py` (Task-2 point checks + identity residual), `step_12.py` (naïve interval artifact, transparency),
`step_13.py` (mean-value mesh certification + derivative validation), `step_14.py` (saddle containment),
`step_15.py` (worst-subcell refinement), `step_16.py` (G-scans + certify_H3 verbatim), `step_17.py`/`step_18.py`
(step_18 timing/logic rerun), `step_19.py` (verbatim step_18 subprocess + json diff), `step_20.py`
(tight-vs-coded + slack decomposition), `step_22.py` (full-digit saddle margins + G sign structure),
`step_01`–`step_10.py` (file reading / line quoting), `step_21.py` (listing), `result.md`.

**Residual risks.** (i) The mean-value mesh certifies containment over the continuum of every subcell, but is
a 64×64 mean-value enclosure rather than a symbolic proof — mitigated by the 256×256 refinement agreeing to
3.1e-8 and by the point checks. (ii) `ψ(t) ≤ (√2/3)t³` (upper box, all t > 0) is evidenced by
`fixed/verify_mS.py` (t ∈ [10⁻¹⁰, 10¹⁰], min margin > 0) and certified on the cover by the mesh; a one-line
symbolic proof would polish it. (iii) D13t's R6 reruns of `step_10_J192.py` and the sweep cells remain
outside this note's scope.
"""
open('result.md','w').write(md)
print("result.md written,", len(md), "chars")
import os
print(os.path.abspath('result.md'))
