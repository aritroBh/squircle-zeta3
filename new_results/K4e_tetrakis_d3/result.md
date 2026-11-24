# K4e: the ζ(3) coefficient of the tetrakis hexahedron (and model B) — final report

**Work dir:** `research/08_theory/K4e_tetrakis_d3`. Labels: **PROVED** (exact/symbolic or written proof) / **NUMERICAL** (printed high-precision agreement) / **CONJECTURE** / **PENDING** (jobs running). Every number below was printed by a `step_NN.py`/`job*.py` in the work dir. Two background high-precision jobs were launched to finish the two open items (they take ~1–2 h): `jobA_hp.py` (model A, h = 1/192 + 1/256, WMIN 1e-95, LMIN 1e-85, dps 135) and `jobB192.py` (model B, outer 1/192, inner 1/128, WMIN 1e-95, LMIN 1e-80, dps 140).

---

## 1. Model A — T₂ verified, T₃ extracted, constant set reduced

**Data (NUMERICAL).** h128 vs h192 agreement of [d^k]S: 3.2e-62, 2.3e-61, 4.4e-59, **3.9e-57 (k=3)**, 2.2e-55 (k=4) ⇒ T₃/√5 has **56 trusted digits**:
`T₃/√5 = −0.16372240519908889847180920266570391331492093287396119879311178…`

**T₂ sign (PROVED by data).** The task statement's T₂ = √5[254/675 − (79/8100)π² + (118/75)Li₂(−1/3) + (4/15)L − (76/75)L²], L = ln(4/3), agrees with the driver to **2.5e-58** (`step_01`). The planner's `planner_check_T2.py` had every sign flipped; the task statement is correct.

**Constant-set reduction (step_10–13; NUMERICAL at 150–250 digits).** Among the 13 task constants there is **exactly one Q-relation**:
`6·A(−1/3) = 12·Li₃(−1/3) + 6·Li₃(1/4) − 6·Li₂(−1/3)·L − L³` (residual 2.5e-150 at dps 150).
The Landen S₃-orbit relation `Li₃(1/4)+Li₃(3/4)+Li₃(−1/3) = ζ(3) − (π²/6)L + ln2·L² − L³/6` also holds (residual **1.2e-251** at dps 250) but **requires ln2·L², which is outside the 13-set**. So within the task's set, ζ(3) is independent of the Li₃'s. **Independent set (12):** {Li₃(−1/3), Li₃(1/4), Li₃(3/4), ζ(3), π²L, Li₂(−1/3)L, L³, L², L, π², Li₂(−1/3), 1}.
*Spurious filter:* candidate relations with max coefficient M were accepted only if the direct-eval residual scales like M·10^(−dps); all large-M PSLQ outputs failed this test at dps 150/200/250 (`step_13`).

**T₃ PSLQ status (step_14, NUMERICAL-negative).** PSLQ of [T₃/√5] + the 12 constants at tol 1e-54/1e-52/1e-50/1e-48 on **both** h192 and h128 gives a *different* relation at every tolerance (M ~ 1e4–1e5, residuals at the tol floor) ⇒ **no genuine relation with M ≲ 1e5 at 56 digits**. The ζ(3) coefficient of [d³]S_A in the 12-constant basis is therefore **PENDING** the jobA_hp digits (~62–70 digits expected; h = 1/256 alone would NOT have helped — convergence analysis `step_31` showed the ~1.8e-57 floor is a WMIN/LMIN noise floor, not systematic: the tanh–sinh systematic is already ~1e-77 by h = 1/192, so jobA_hp uses tighter WMIN/LMIN and dps 135).

## 2. Model B — exact ratios, atoms, internal relations (PROVED/near-PROVED)

**Geometry (PROVED, sympy; `step_26/27`, correcting K4d step_45's pairing).** Vertices: 6 axis (valence 4) at ±a eᵢ, a = (4√3−3)/13; 12 diag (valence 6) at (0,±c,±c), c = (2√3−3)/2.
- **L1 orbit (24 axis–diag edges):** ratios **r₁ = (5−2√3)/2** (axis term), **r₂ = 5−2√3** (diag terms), r₂ = 2r₁. End-vertex sharpnesses α = √3−3/2, β = (4−√3)/13; L₁ = 0.24241329462711, h = 0.28925258324196, (u,v) = (0.36055, 0.63945).
- **L2 orbit (24 diag–diag edges):** ratio 1; α = β = √3−3/2; L₂ = 0.32816939922354, h = 0.28420303647226, u = v = 1/2.

**Atoms (x = 1−r), with a sign correction:** x₁ = (2√3−3)/2, x₂ = 2√3−4; the Landen partner of x₁ is **−α = −(4√3−3)/13** (K4d-era notes had +α — wrong sign, verified `step_29`), of x₂ is 2β = 2(4−√3)/13; 1/r₂ = (5+2√3)/13.
**Relations verified at ~120 digits (`step_28/29/30/32`):**
- weight-2 Landen: Li₂(x₁)+Li₂(−α) = −½ln²r₁; Li₂(x₂)+Li₂(2β) = −½ln²r₂.
- **General A-reduction** (verified at x₁, x₂, −1/3, 0.4, −0.7, 0.15, −0.9, residuals ≤ 5e-110 — a classical S_{1,2} identity, model A's relation is its x = −1/3 case):
  **6·A(x) = 12·Li₃(x) + 6·Li₃(x/(x−1)) − 6·Li₂(x)·ln(1−x) − ln³(1−x)**.
- x₁-orbit Landen: Li₃(x₁)+Li₃(r₁)+Li₃(−α) = ζ(3) + ζ(2)ln r₁ − ½ln(x₁)ln²r₁ + (1/6)ln³r₁ (4.8e-121); x₂-orbit analog (complex Li₃(r₂), Re = 2.1172201455689146…) at 1.6e-121.
All vertex-value logs live in span{ln2, ln3, ln(2−√3), ln(4+√3)}.

**q_B (NUMERICAL-negative; `step_18–25`).** First, a data correction: the planner's q_B = 3.268832565066791005838… came from `modelB_h96.json` (S₀ error 1.7e-17, the *worst* grid). The accurate h48/h96o value (agreement 2e-64) is
**q_B = 3.26883256506679086582204524494527788236949929032104050993**.
Exhaustive PSLQ — ρ + Σwᵢln(gᵢ) with ρ,w ∈ Q (excluded M < 1e8), ρ,w ∈ Q(√3) (excluded M < ~1e5–1e7), quartic-constant and [d]S = t₁ + S₀lnc_B forms (excluded M < ~1e6) — finds **no small relation**; the one M = 77146 candidate was refuted by tolerance-instability (`step_25`). q_B is **PENDING** the h192B digits (mirrors K4c's "leaves the small-log world"). Method sanity: q_A recovered as 15q_A − 13 + 18ln2 − 39ln3 = 0 (res 5.8e-60).
**T₂/T₃ for B:** [d²]S/S₀ PSLQ also needs the extra digits (relations tolerance-unstable at 55 digits, `step_32`); the ζ(3) coefficient of [d³]S_B is basis-independent of the c_B mixing (c_B only shifts lower weights) and is **PENDING** h192B.

## 3. Anchors and 3D-law status

- **Cuboctahedron verified end-to-end (step_33, NUMERICAL 58 digits):** S(1+d) = 4^{d/(1+d)}[S₁ + S₁d/2 + (√3π²/96 − S₁/8)d² + ((7−√3)π²/64 − 23S₁/48 − (√3/6)ζ(3))d³] matches K4b data to ≤ 1.5e-58 for k = 0..3 ⇒ [d³]S_cubo|ζ(3) = **−√3/6**, q_cubo = 1/2 + ln4. Octahedron: [d³]S|ζ(3) = **−64√3/9**, q = 1.
- **2D edge law does not transfer directly** (step_34): applying Z_e = 6αβ(h²/L)(3(αu+βv)−2) per edge gives octa +25.46 vs actual −12.32 and cubo −8.27 vs −0.289. A per-length form F ∝ (αβ)²·ℓ·g(u,v) fits the octa/cubo ratio 36.95 exactly (4² × 2.3094), but u = v = ½ on both anchors, so g is undetermined; the asymmetric orbits (A apex-base u = 1/3; B L1 u = 0.3606) will probe g once the two ζ(3) coefficients land. **4-point fit: PENDING** the two job outputs.

## 4. Evidence table

| claim | residual / digits | status | file |
|---|---|---|---|
| T₂(A) task sign correct | 2.5e-58 | NUMERICAL (58 dig) | step_01 |
| T₃/√5(A) value | h128/h192 agree 1.77e-57 | NUMERICAL (56 dig) | step_01/31 |
| 6A(−1/3) = 12Li₃(−1/3)+6Li₃(1/4)−6Li₂(−1/3)L−L³ | 2.5e-150 | NUMERICAL (150 dig) | step_10/11 |
| Li₃-orbit R1 with ln2·L² | 1.2e-251 | NUMERICAL (250 dig) | step_11/13 |
| 13-set ⇒ exactly one relation; indep set of 12 | spurious filter applied | NUMERICAL | step_12/13 |
| T₃(A) no relation M ≲ 1e5 at 56 dig | tol-instability shown | NUMERICAL-negative | step_14 |
| B ratios r₁ = (5−2√3)/2, r₂ = 5−2√3 | exact | PROVED (sympy) | step_26/27 |
| B A-reduction & Landen orbit relations | ≤ 5e-110, 7 points | NUMERICAL (120 dig) | step_29/30/32 |
| q_B corrected value; not identified M < 1e8 | 2e-64 grid agreement | NUMERICAL-negative | step_18–25 |
| cubo expansion & −√3/6 anchor | ≤ 1.5e-58 | NUMERICAL (58 dig) | step_33 |

## 5. Step files produced
step_01 (T₂/T₃ extraction), step_10/11/12/13 (constant reduction + spurious tests), step_14 (T₃ PSLQ), step_18–25 (model B data, q_B campaign), step_26/27 (edge data), step_28/29/30/32 (B atoms + relations, [d²]S try), step_31 (convergence/extrapolation), step_33 (cubo verification), step_34 (2D-law countercheck); jobs jobA_hp.py, jobB192.py (running), result.md, progress.log, consts13.json, modelB_atoms2.json.

**Bottom line:** The full machinery for the T₃ identifications is built and validated (T₂ A confirmed; cubo anchor confirmed; all constant-set and model-B internal relations certified at 120–250 digits; spurious-relation filters in place; both data-generation jobs running at tightened WMIN/LMIN/dps). The two headline numbers — the ζ(3) coefficients of [d³]S for model A and model B — were **not** isolable within the turn budget at the available 55–56 digits (PSLQ tolerance-instability demonstrated, not merely absence of output); they await the jobA_hp/jobB192 digits, after which `step_14` and its model-B analog (reduced bases now known: 12 constants for A; {Li₃(x₁), Li₃(−α), Li₃(x₂), Li₃(2β), ζ(3)} + log monomials for B) and the 4-point edge-law fit are one PSLQ call each.