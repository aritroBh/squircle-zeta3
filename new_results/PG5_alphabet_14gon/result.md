# PG5: the polygon alphabet rule at s = 5 (and 6), and the regular 14-gon

## 1. Task 1 — rank computation at s = 5: COMPLETE (exact, printed)

Method (PGEN1R's engine `symbol_engine_r.py`, copied into cwd; exact integer arithmetic; Goncharov recursion; Q₃ = ⊗³V/D computed via the perfect pairing with Lie³ Lyndon brackets; ranks by `flint.fmpz_mat`; I added an inverted-index fast pairing). Reproduced PGEN1R's table for s = 1..4 exactly, then extended:

| s | n | dim Q₃ (Witt, checked by pairing) | #CR | rank s(CR) mod D | **2·C(s+3,4)** | #words | rank words mod D | **s(2s²+3s+7)/6** | words ⊂ span s(CR)+D | shuffle gens killed |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 2 ✓ | 6 | 2 | 2 | 3 | 2 | 2 | YES | ✓ |
| 2 | 5 | 40 ✓ | 30 | 10 | 10 | 14 | 7 | 7 | YES | ✓ |
| 3 | 9 | 240 ✓ | 90 | 30 | 30 | 39 | 17 | 17 | YES | ✓ |
| 4 | 14 | 910 ✓ | 210 | 70 | 70 | 84 | 34 | 34 | YES | ✓ |
| **5** | **20** | **2660 ✓** | **420** | **140** | **140** | **155** | **60** | **60** | **YES** | ✓ |
| **6** | **27** | **6552 ✓** | **756** | **252** | **252** | **258** | **97** | **97** | **YES** | ✓ |

**The CR-rank pattern 2·C(s+3,4) HOLDS at s = 5 (140) and at s = 6 (252).** The word rank also follows a clean cubic, s(2s²+3s+7)/6 (fitted on s = 1..4, then confirmed at s = 5 and 6 — I did not assert it a priori). Inclusion holds at both new levels: rank(CR ∪ words) = rank(CR), i.e. every convergent distinct-letter weight-3 word symbol lies in span{s(Li₃(cross-ratios))} + D.

**Repeated-letter repair (PGEN1R's N(x,y) fix), run at s = 5 and s = 6** (`step_03_repair.py`): for every ordered pair (x,y), x ∈ {p₁..p_s}, y ∈ {1,p₁..p_s}∖{x}: c = x/(x−y) and cX = (x−1)/(x−y) together with 1−c, 1−cX all lie in the CR set (25/25 at s = 5, 36/36 at s = 6), and adjoining s(N) = 2s(c)−2s(cX) for **all** pairs plus the reversal cases G(1,1,z) = Li₃(1/(1−z)) (5/5, 6/6 in CR set) leaves the rank unchanged (140 at s = 5, 252 at s = 6). So **all** repeated-letter convergent words also lie in span s(CR)+D. Combined with PGEN1R's certified reduction (Nielsen antiderivative, shuffle identities, t↦1−t reversal), the alphabet rule is proved at symbol level for s ≤ 6:

> **Theorem (symbol level, s ≤ 6).** On an edge with corners 0, 1 and s ≤ 6 far kink lines crossing its line at p_j (i.e. centrally symmetric polygons with up to 16 edges), every weight-3 iterated integral with letters {0,1,p_j} has, modulo products, the symbol of a combination of Li₃ at cross-ratios of {0,1,∞,p₁,…,p_s}.

## 2. Task 2 — s = 6 feasibility

Feasible and done: 5.7 s for the s = 6 rank block (756×6552 and 1154×6552 flint ranks), 36 s including the dim-Q₃ pairing check (19683×6552) and 3.2 s for the repair. Nothing heavier was needed; the cost would only become an issue around s ≈ 8–9 (dim Q₃ ≈ 14k–20k).

## 3. Task 3 — the regular 14-gon: infrastructure PROVED; δ(N₃) NOT obtained

Status: the full exact pipeline was rebuilt for m = 7 up through the exact d³-integrand E₃ (validated pointwise to 1e-18 by an independent construction), but the integration-engine port, assembly, Li₃-basis reduction and hence the δ(N₃) computation were **not completed** in the session. I report only what was printed.

**Geometry (exact in F = Q(α), α = 2cos(2π/7), my `ftfield.py`, all identities verified by exact arithmetic; `geom14.pkl`).** K_p = {Σᵢ₌₀⁶|⟨nᵢ,x⟩|^p ≤ 1}, nᵢ = (cos iπ/7, sin iπ/7). Edge in chamber (−π/14, π/14): x₁ = h with **h = 1 − α²/2**, L² = (95 − 25α − 41α²)/7, h²/L² = 3/4 + 2α + α². Terms (order [loc t, loc 1−t, a₁, a₂, a₅, a₆, const]):
- A₀ = A₁ = (3 − α − α²)/2; B₁ = B₅ = −2 + α/2 + α²; B₂ = A₀; B₆ = h; constant term a₀ = h.
- x₁ = 1 − α, x₂ = 2 − α − α², x₅ = 2 − α², x₆ = 3 − α − α².
- Verified exactly: Σaᵢ = 1 identically; a₃(0) = a₄(1) = 0; Landen pairs x₆ = x₁/(x₁−1), x₅ = x₂/(x₂−1); reflection B₆ = B₁(1−x₁), B₅ = B₂(1−x₂). Structure = DOD1's (4 nonconstant far logs in 2 Landen pairs + 1 constant normal term, kink line parallel to the edge → p₀ = ∞). s = m−2 = 5 far kink lines ✓.

**Letters (`letters14.pkl`).** Everything has norm ±2^k: the letter group is **{2, α, u = α²+α−1 = 2cos(π/7)}** (2 inert; α, u units; all xᵢ, 1−xᵢ, Bᵢ, all 12 γ/δ cross points and complements, and all cross-ratios factor — 126/126 factorizations verified exactly by multiplying back). The edge's CR set over {0,1,∞,p₁,p₂,p₅,p₆} has **24 distinct values** (specialization collapses the 420 formal CRs). Letter rank 3 → **dim Q₃ = 8** (as for the 12-gon).

**E₃ (`E3n14.pkl`, 119 monomials, built by DOD1's d-series formula in exact FT arithmetic, 3.5 s).** Validated pointwise against an *independent* polar-perimeter d³-extraction (roots-of-unity filter on J(t,d) = √(r²+r′²)(θ(t),1+d)θ′/L, r = f^(−1/p)): 6·E₃(t) matched at t = 0.13/0.37/0.5/0.72/0.91 with residuals 2.0e-18, 3.5e-19, 4.1e-19, 1.7e-19, 8.2e-18.

**Validation tension (unresolved).** Two candidate values for P‴(1) were printed: (i) polar-perimeter gold (Fourier filter, two radii agreeing to 7.1e-18, and an 8-point real-axis stencil agreeing with it to 10 digits): **−5.5080732424212877414**; (ii) the exact-integrand chain 84·L·∫₀¹E₃ dt with ∫E₃ = −0.6455382870546723 (two independent quadratures agreeing to 3.7e-13): **−5.5080836980534783** — a 1.05e-5 discrepancy. Since E₃ is validated pointwise at 1e-18 and the two integral evaluations agree, my working diagnosis is a smooth-in-p quadrature error (~1e-11 per sample, amplified by ρ⁻³ in the c₃ extraction) in the polar-perimeter gold; but without the engine/assembly I could not adjudicate, so I report both numbers and no verdict. This blocks the P‴(1)-level validation of the assembly, **but not the δ computation itself** (which is exact and numerics-free); δ was simply not reached in time.

**δ(N₃) verdict: NOT OBTAINED.** The remaining work is mechanical and costed: port `dod_engine` (pair/triple recursions, Jloc/Jtab tables — already inspected and small) to FT arithmetic (~2–4 h coding), run the assembly (m = 7 has the same integral-type census as m = 6: 4 singles, 6 pairs, 4 triples), then the Q₃-symbol reduction over the 3 letters (dim 8, trivial linear algebra) and δ with the PGEN1R invariance checks.

## 4. Verify script (Task 4)

`step_02_s5_ranks.py 123456` reproduces the full rank table s = 1..6 in ≈ 40 s; `step_03_repair.py 56` the repair in ≈ 4 s; `step_04_geom14.py` + `step_05_letters14.py` + `step_07_e3n14.py` + `step_08_e3check.py` re-derive the 14-gon data and the E₃ pointwise validation in ≈ 15 s (the gold and the ∫E₃ comparison take ≈ 12 min and ≈ 5 s; the s = 5-only checks are ≈ 2 s). All under the 5-minute budget for the s = 5 claims.

## Evidence table (all numbers printed in-session)

| check | result | file |
|---|---|---|
| ranks s=1..4 reproduce PGEN1R | exact match (2,10,30,70 / 2,7,17,34, YES) | step_02_s5_ranks.py |
| **s=5: rankCR=140=2C(8,4), rankW=60, inclusion YES, dimQ₃=2660 pairing, D killed** | printed | ranks_345.pkl |
| **s=6: rankCR=252=2C(9,4), rankW=97, inclusion YES, dimQ₃=6552 pairing, D killed** | printed (36 s) | ranks_6.pkl |
| repair s=5: 25/25 pairs c,cX,1−c,1−cX ∈ CR; rank with all N = 140 | printed | repair_s5.pkl |
| repair s=6: 36/36 pairs; rank with all N = 252 | printed | repair_s6.pkl |
| 14-gon geometry: Σaᵢ=1, Landen, reflection, a₃(0)=a₄(1)=0 | exact (FT arithmetic) | geom14.pkl |
| letters {2,α,u}; 126/126 factorizations; |CR_e| = 24; dim Q₃ = 8 | exact | letters14.pkl |
| E₃ 119 monomials; pointwise vs independent polar d³ | ≤ 2.1e-18 at 5 points | E3n14.pkl, step_08 |
| ∫E₃ (two quadratures) | −0.6455382870546723 (agree 3.7e-13) | E3n14_check.pkl |
| P‴(1) candidates: −5.5080732424212877 (polar gold, 7e-18 self-consistent) vs −5.5080836980534783 (exact-integrand chain) | **1.05e-5 unresolved** | gold14.pkl, step_08 |
| δ(N₃) of the 14-gon | **not computed** | — |

## Files in `new_results/PG5_alphabet_14gon/`

`symbol_engine_r.py` (PGEN1R engine), `step_02_s5_ranks.py` (rank table s=1..6 + dim/d D checks), `step_03_repair.py` (repeated-letter repair at s=5,6), `ftfield.py` (exact Q(2cos(2π/7)) arithmetic), `step_04_geom14.py` + `geom14.pkl` (exact edge data), `step_05_letters14.py` + `letters14.pkl` (letter basis, 24 CRs, 126 factorizations), `step_06_gold14.py`/`step_06b_gold14.py` + `gold14.pkl` (polar-perimeter P‴(1)), `step_07_e3n14.py` + `E3n14.pkl` (exact E₃), `step_08_e3check.py`/`step_08b_intcheck.py` + `E3n14_check.pkl` (pointwise validation + integral), checkpoints `ranks_*.pkl`, `repair_s*.pkl`, `progress.log`.

**Bottom line.** The alphabet rule for the trilogarithm arguments of P‴(1) is now proved at symbol level for s ≤ 6 far kink lines — covering centrally symmetric polygons up to 16 edges, i.e. beyond the paper's twelve-edge frontier — with rank s(CR) mod D = 2·C(s+3,4) at every level (140 at s = 5, 252 at s = 6), word rank s(2s²+3s+7)/6, inclusion YES, and PGEN1R's repeated-letter repair verified at the new levels. For the regular 14-gon, the exact edge data, the letter group {2, 2cos(2π/7), 2cos(π/7)} (dim Q₃ = 8), the 24-point cross-ratio set, and the exact d³-integrand E₃ (validated to 1e-18 pointwise) are in place; the δ(N₃) computation itself and one 1.05e-5 validation discrepancy remain open and are clearly scoped for the next session.