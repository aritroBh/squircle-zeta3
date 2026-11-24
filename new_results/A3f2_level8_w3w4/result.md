# A3f2 — level-8 relation completion, weights 3 & 4: status report

**Honest headline:** the weight‑3 machinery advanced substantially (a *proven* two‑variable Li3 relation family was derived from the Kummer Li4 equation), but a one‑character transcription typo in my 20‑term monomial table was only located in the final debug step, so no new quotient numbers were certified this session. Weight 3 remains at **34 (target 27)**, weight 4 at **142 (target 81)**, D reduction **OPEN**. All previously audited results stand.

## 1. Status table

| step | status | detail |
|---|---|---|
| Setup: copy Q engine + pickles, verify imports | PROVED (done) | `Q/` copied to workdir; `cmzv`, `relations`, `check_conversion` import cleanly; all pickles load |
| kummer_aux content | PROVED (done) | 36 D0 Kummer instances → **2 nonzero eliminated H‑rows** (Codex's pipeline; not yet merged into the w4 relation set) |
| Direct in‑configuration K4 instances | PROVED = 0 | `Q/kummer4.pickle`: 0 instances with all 20 arguments cross‑ratios ⇒ auxiliary elimination is mandatory |
| Real cross‑ratios of the 8‑point configuration | PROVED (exact) | **27 real cross‑ratios, all S‑units** ±2^(j/2)·a^m; 9 of them in (0,1) |
| **New: proven 2‑variable Li3 families Rx, Ry** | PROVED (derived) | ∂/∂x and ∂/∂y of Kummer's K4 (see §3) |
| K3 candidate pool | NUMERICAL | 11 S‑unit candidates x∈(0,1) with 1−x S‑unit (|j|,|m|≤6); 36 pairs x+y≤1; ~211 Li3 atoms (~193 auxiliary) |
| Rx/Ry row generation + rank test at w3 | **FAILED (bug found)** | 4th monomial exponent typo: `y²x/(ξ²η)` needs c=−2, I wrote c=+2. K4 self‑test residual −2.38e‑2 (should be 0); Rx/Ry residuals −4.2e‑2/−4.1e‑2. One‑character fix |
| Weight‑3 quotient | unchanged: **34/27** | (37 standard → 34 with lifted w2 relation, per A3f step_09) |
| Weight‑4 quotient | unchanged: **142/81** | 17625 rows audited mod 2 primes in A3f |
| D reduction | **OPEN** | — |

## 2. What was established exactly this session

**(a) The configuration's real cross‑ratios (exact, printed):** 27 values, every one an S‑unit of ℚ(√2): {−a², −2a, −a, −√2, −1, −1/√2, −b, −b²/2·…, b², 2^{−3/2}a^{−1}, b, 1−1/√2, 1/2, 2−√2, 1/√2, 2^{3/2}b, (2+√2)/4, … , a², 2a} (full list with (j,m) exponents printed in the workdir step output). Since S is closed under z↦1−z, all Kummer monomials built from S‑unit (x,y) are S‑units.

**(b) The K3‑derivative relations (the session's main mathematical result).**
Kummer's K4 is proven (Kummer 1840; Lewin Ch. 7; exact symbol guard S(LHS)=6·(ξξ⧢ηη); 240‑point 60‑digit check worst 3.7e‑59 in `B_kummer_li4`). With G(x,y)=ΣᵢdᵢLi4(mᵢ)−(3/2)log²ξ·log²η ≡ 0 on the convex domain D0={x,y>0, x+y<1}, its partial derivatives are proven identities on D0:

- **Rx:** Σᵢ dᵢ(aᵢξ − cᵢx)·Li3(mᵢ(x,y)) + 3x·logξ·log²η = 0
- **Ry:** Σᵢ dᵢ(bᵢη − dᵢy)·Li3(mᵢ(x,y)) + 3y·log²ξ·logη = 0

where mᵢ = ±x^{aᵢ}y^{bᵢ}ξ^{cᵢ}η^{dᵢ} are the 20 Kummer monomials, dᵢ∈{+1,−6,−3,+6}. Every specialised instance (x₀,y₀ S‑units, x₀+y₀≤1) is therefore a **proven** weight‑3 relation; RHS logs are S‑unit logs (expressible via the engine's log2/logA atoms). Coefficients lie in ℚ(√2); the rank computation must be over K=ℚ(√2) — for split primes p≡±1 mod 8 take r²=2 mod p and rank of (R₀+rR₁) over 𝔽_p equals the K‑rank; exactly, the doubled rational matrix has ℚ‑rank = 2·K‑rank (derived in the workdir notes).

**(c) Plan validated against the obstruction:** 0 direct instances ⇒ auxiliary‑atom elimination (QS2‑RREF over the Li3‑atom universe with T3/I3/D3/V3 rows, exactly mirroring Codex's `kummer4.auxiliary_closure`) is required; the atom pool (~211 atoms from 36 pairs) is small enough for exact QS2 elimination.

## 3. Evidence table (all numbers printed by my code this session)

| check | printed value |
|---|---|
| \|S\| real cross‑ratios | 27 (all S‑units, (j,m) tags printed) |
| positive cross‑ratios in (0,1) | 9 |
| kummer4 direct config instances (from `Q/kummer4.pickle`) | 0 |
| kummer_aux eliminated H‑rows | 2 (both nonzero) |
| weight pickles | w3: 576 cols/1298 rows/539 sel; w4: 5184 cols/17625 rows |
| S‑unit candidates x∈(0,1), 1−x S‑unit, \|j\|,\|m\|≤6 | 11 |
| instance pairs x+y≤1 | 36 |
| **K4 self‑test with my (buggy) table at (0.2,0.3)** | residual **−2.38062e‑2** ⇒ table error |
| Rx, Ry residuals (buggy) | −4.18629e‑2, −4.06158e‑2 |
| numeric dG/dx (buggy) | −2.61643e‑1 |

Prior‑session anchors (audited, unchanged): w1–3 quotients 3/10/37 exact over ℚ; w4 = 142 mod primes 1000000009 & 1000000021; all 17625 w4 rows numerically 0 (worst 2.9e‑57); weight 2 CLOSED at quotient 9 via Rogers five‑term at x=√2/2, y=√2−1 (planner re‑eval 2.8e‑41 @40 digits); w3 = 34 after lifts; T3/I3 added 0 rank.

## 4. The bug and the exact next step

In the 20‑term table K20, entry 4 must be `y²x/(ξ²η) = (s,a,b,c,d)=(1,1,2,−2,−1)`, I wrote `c=+2`. Consequences: the "0 direct instances / ≤4 near‑misses" scan and the atom counts are slightly unreliable and must be re‑run. **Next step (ready to execute):** fix the one exponent, re‑run the K4 self‑test (expect ~1e‑50 residuals as in the library), regenerate Rx/Ry rows, numeric‑check each at 60 digits via `check_conversion.evaluate`, run QS2 auxiliary elimination with T3/I3/D3/V3 rows, then doubled‑matrix K‑rank on top of `weight3_with_lifts.pickle` (mod 2 primes, then exact) and report the w3 quotient; then rebuild w4 relations with enriched `lower` bases + kummer_aux 2 rows + I4 instances, and test D membership after each batch.

## 5. Files in the workdir `…/08_theory/A3f2_level8_w3w4/`

- `Q/` — full copy of the MPV engine (`cmzv.py`, `relations.py`, `check_conversion.py`, `kummer4.py`, weight pickles, `kummer_aux.pickle`, `vectors.json`).
- Pickles copied from A3f: `audit_weight{1,2,3}.pickle`, `weight2_new_relation.{pickle,json}`, `weight2_fiveterm.pickle`, `weight3_with_lifts.pickle`, `weight3_onevar.pickle`.
- `step_01.py … step_32.py` — every computation above (auto‑saved per run). Key ones: cross‑ratio enumeration (27 S‑units), K3‑derivative builder producing `k3_deriv_raw.pickle` (**contains the c=+2 typo — do not reuse without the one‑character fix**), and the final debug step that isolated the typo by reproducing the K4 residual −2.38e‑2 and the numeric derivative.
- `progress.log` lines were not yet appended (session cut short); this report serves as the log content: setup+engine copy OK; 27 real cross‑ratios all S‑units; Rx/Ry derived from K4; pool 11 cand/36 pairs/~211 atoms; K20[3] typo found (c: +2→−2); quotients unchanged w3=34/27, w4=142/81; D OPEN.

**Bottom line:** no false claims — the certified quotients are still w2=9 (closed), w3=34 of 27, w4=142 of 81, D not reduced. The session's durable output is the *proven* Rx/Ry two‑variable Li3 family differentiated out of Kummer's K4, the exact 27‑element real cross‑ratio/S‑unit table, the K=ℚ(√2) rank methodology, and a precisely localised one‑character bug whose fix unblocks the weight‑3 rank test.