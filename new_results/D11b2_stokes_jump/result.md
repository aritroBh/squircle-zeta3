# D11b2 — Stokes-jump check at the diamond end of the squircle

Working directory: `research/08_theory/D11b2_stokes_jump/`. Data from D11b: `step2_num450_N180.json` (f₀…f₁₈₀ at 450 digits), `pade1_9090.json`. All numbers below were printed by code run in this directory (step files listed at the end). Three heavy integration runs hit the 1500 s cap (stdout lost on kill); the consequences are stated honestly in §4.

---

## 1. Where are the Borel singularities? The "poles below π/2" mystery — resolved

### 1.1 Direct test of the Borel coefficients (NUMERICAL, confirming D11's PROVED asymptotics)

With b_n = f_n/n! from the 450-digit data, D11's proved formula gives
b_n ~ pat(n)(2/π)ⁿ n^(−3/2), pat(n) = √2 cos(π(3−2n)/4) = [−1,1,1,−1] (|pat|≡1).

**Root test** (`step_08.py`): |b_n|^(−1/n) decreases toward π/2 with the predicted slow log-correction:

| n | \|b_n\|^(−1/n) | /(π/2)−1 |
|---|---|---|
| 60 | 1.74052204504547083 | 0.1081 |
| 100 | 1.6832811494478092 | 0.07161 |
| 140 | 1.65627554492266743 | 0.05442 |
| 180 | 1.64030610749119712 | 0.04425 |

Fit |b_n|^(−1/n)/(π/2) = 1 + a·ln n/n + c/n on n∈[120,180]: a = 1.37416 (theory 3/2), c = 0.82833, rms 2.96e−6. ✓

**Domb–Sykes, period 4** (pat(n+4)=pat(n) cancels the oscillation; r_n = (b_n/b_{n+4})^{1/4} = R(1+a₁/n+…)): 6-point fit ns=[120,…,170], p=4:

> R = 1.570796327847109026945 vs π/2 = 1.5707963267948966, **error 1.05e−9 (9 digits)**;
> a₁ = 1.4999995 (theory 3/2, 7 digits) ⟹ coefficient power n^(−3/2) ⟹ **square-root (exponent 1/2) tips**.

The sign pattern pat(n) = Re[e^{3iπ/4}(−i)ⁿ·2/√2·…] fixes the tip locations at exactly **t = ±iπ/2** (terms (∓i)ⁿ(2/π)ⁿ ↔ Li_{3/2}(±2it/π), branch at ∓2it/π=1). Given D11's proved asymptotic, radius = π/2 with the only boundary singularities at ±iπ/2 is **PROVED**; the 9-digit Domb–Sykes agreement is the independent NUMERICAL confirmation.

### 1.2 The poles at ±1.4516i/±1.4530i are a root-finding artifact — the true Padé poles all lie beyond π/2 (NUMERICAL)

* The Padé denominators have coefficients spanning ~10⁹⁵ (b_n ranges over 95 orders of magnitude). Roots extracted by `numpy.roots` (double precision) — as used in D11b — are meaningless in the inner region. I reproduced the fake values: PB1[90/90] "poles" at −0.05429+1.33844i (|t|=1.33954), +0.04709+1.34884i (1.34966)…; PB1[70/70] at −0.00353+1.44383i (1.44384). Even den(root)≈5e−11 "verification" is vacuous at this coefficient scale. The builds themselves are stable (max |den coefficient difference| between dps-350 and the saved dps-400 builds = 1.23e−55).
* **Newton refinement at dps=100 moves every one of the 46 inner candidates to |t| ≥ 1.57189 > π/2** (`step_10/11`): innermost true pole of PB1[90/90] at t = 0.008639114665382 + 1.571867252489j, |t| = **1.57189099296**; then 1.57213, 1.57216, 1.57220, … marching up the cut, straddling the axis (Re ≈ ±0.01) — textbook cut discretization. PB1[70/70] (refined, dps=80): innermost |t| = **1.57226548358**. As N grows 70→90 the innermost pole moves 1.57227 → 1.57189, i.e. **converges down toward π/2 = 1.57079633 from above** — exactly the expected behavior for a branch-point cut.
* Calibration on the 1-cut model Li_{3/2}(2it/π) (whose only singularity is the tip at −iπ/2): numpy likewise reports fake poles at |t| = 0.60–0.80, while properly refined poles are all on the cut ray at |t| = 2.3456, 2.5709, 2.8526, 3.2086 (z = 2it/π = 1.4932, 1.6367, 1.8160, 2.0426 > 1), as Stieltjes theory demands. The 2-cut model with the exact D11 leading coefficients (pat phases) reproduces the same fake-inner-pole pattern as the real data (1.32–1.45), proving the phenomenon is generic to this Padé+root-finding combination, not to the data.

**Answers to the task's questions, plainly:**
* The true singularity **is** at ±iπ/2 (PROVED via D11; confirmed to 9 digits by Domb–Sykes). There is **no** singularity at ±1.4516i.
* The reported poles are **not** Froissart doublets (their apparent residues were ~0.015, not tiny, and they are stable under rebuild precision at fixed N) and **not** an ηs-vs-t scaling artifact (both PB and PB1 show them; the Borel variable is the correct one). They are **spurious roots produced by double-precision polynomial root-finding on an extremely ill-conditioned Padé denominator**; honest refinement places all poles on the cuts beyond ±iπ/2.

## 2. The Stokes jump

### 2.1 Prediction and normalization (PROVED, given D11's asymptotic)

From B(t) ⊃ −√(2π)e^{3iπ/4}(1+2it/π)^{1/2} at t = iπ/2 and the Hankel contour, with J = S_left − S_right (left = ray arg s = +δ, i.e. Re t < 0 side of the upward cut) for the Borel sum of **F itself**, F = √2/(1+η) + (1/η)∫e^{−t/η}B₁(t)dt:

> **J(ε) = 2√π e^{−3iπ/4} ε^{1/2} e^{−π/(2ε)} (1 + O(ε))**, |C| = 2√π = 3.544907702, **b = 1/2**, phase **e^{−3iπ/4}**.

Normalization check: the task's stated "b = 3/2, phase e^{+3iπ/4}" is **not** the F-normalization. The ∫e^{−y/ε}y^{1/2}dy ∝ ε^{3/2} from the square-root tip is multiplied by the Borel 1/η = 1/(iε), giving **b = 1/2 for F** (b = 3/2 is the jump of the t-integral ηF). The phase e^{−3iπ/4} holds for J = S_left − S_right as defined above; flipping the order conjugates/negates it.

### 2.2 Method and pipeline validation on the exact model (NUMERICAL)

Conformal map w = (1−√(1−u²))/u, u = −2it/π (inverse t = iπw/(1+w²)) unfolds the square-root: 1+2it(w)/π = (1−w)²/(1+w²), so B̃(w) = B(t(w)) is **analytic at w=1**; the two cut sides become the unit arcs w = e^{±iφ}, and
J = (1/(iε)) ∫₀^{π/2} e^{−π/(2εcosφ)} [G(φ)−G(−φ)] dφ, G(φ) = PB̃(e^{iφ}) t′(e^{iφ}) i e^{iφ}.
Built d-series (d_k: |d_k|^{1/k} = 1.051, 1.0468, 1.0434, 1.0406 at k=90,120,150,180 → radius 1), w-Padé [85/85] and [90/90] (poles nearest w=1 at distance ≈0.93; the unfolded tip is clean). mp.quad exhibited a pathology at the w=i endpoint (blew up after 862 s); replaced by 64–96-point Gauss–Legendre.

Model check (exact polylog model with D11's leading coefficients, exact Hankel integral via disc Li_{3/2}(X) = 2πi(ln X)^{1/2}/Γ(3/2)) at ε = 0.2:

> J_model^exact = −4.17208512412e−4·(1+i) — phase exactly **e^{−3iπ/4}** ✓;
> |J| = 5.9001e−4 vs predicted 2√π·√0.2·e^{−7.85398} = 6.1393e−4: ratio **0.961**, consistent with 1+O(ε) ✓.

### 2.3 Real-data values (INCOMPLETE — stated honestly)

The real-data J_arc table (w-Padé [85]/[90] × 5 ε values, plus lateral-sum cross-checks) did not finish: three attempts hit the 1500 s cap (mp.quad pathology; then two full rebuilds at dps 80–120). Progress logs may exist in `step17_progress.log`/`step_17_J_partial.json` but could not be recovered before the deadline. The fit of C, b, phase on the real function therefore remains open; the validated pipeline (`step_16/17`) plus saved `step_12_wpade.json` make it a ~10-minute rerun at dps≤60.

---

## Evidence summary

| Claim | Status | Evidence |
|---|---|---|
| Radius of B is exactly π/2; square-root tips at ±iπ/2 only | PROVED (D11) + NUMERICAL (9 digits) | Domb–Sykes R err 1.05e−9, exponent 1.4999995 |
| "Poles at ±1.4516i below π/2" | NUMERICAL: artifact of double-precision root-finding | Refinement: all 46 candidates → \|t\|≥1.57189 (90/90), 1.57227 (70/70); same fakes in 1-cut and 2-cut models |
| Froissart doublets? Scaling artifact? | NUMERICAL: neither | residues O(0.015); stable vs build dps (1.2e−55); unsubtracted PB also shows fakes |
| J(ε) = 2√π e^{−3iπ/4} ε^{1/2} e^{−π/(2ε)} for F (b=1/2, not 3/2; phase e^{−3iπ/4}) | PROVED from D11's asymptotics (Hankel argument); NUMERICAL on exact model (phase exact, \|J\|/pred = 0.961 @ ε=0.2) | §2.1, §2.2 |
| Real-data fit of C, b, phase | INCOMPLETE (budget) | §2.3 |

## Files produced here
`step_08.py` (root test + Domb–Sykes), `step_09/10/11.py` (Padé pole forensics, refinement; `step_09_pade70/90.json`), `step_12.py` (conformal d-series + w-Padé [85]/[90]; `step_12_wpade.json`), `step_13–15.py` (model calibration, exact model J), `step_16/17.py` (Gauss–Legendre arc-integral pipeline for J; `step_16_J.json`, `step17_progress.log` partial).