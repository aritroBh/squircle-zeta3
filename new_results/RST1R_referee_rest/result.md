# RST1R referee report — sharp rest exponent √(2π) in the squircle decay law

**Under review:** `new_results/RST1_rest_exponent/{proof.md, result.md, verify_rst1.py, rst1_common.py, rst1_common2.py, step_15.py, step_17.py, step_18b.py}`. Claim: for every N = k+1 ≥ 100, |M_k − M_k^main| ≤ 1030·N^{−1/2}·e^{−√(2πN)} (C_raw = 1025.8898389…, arb); corollary C_3 = 1,000,668; item 1 (numerical). Every number below was printed by my own `step_NN.py` files in `new_results/RST1R_referee_rest/` (verify rerun log: `verify_rst1_rerun.log`).

## OVERALL VERDICT: **PROVED AFTER REPAIRS**

The theorem is correct; every load-bearing lemma re-derives independently; the constant C_raw matches to all printed digits; verify_rst1.py reran clean (PASS, 240 s). Two repairs required — both far below the certified margin of 4.11, so the theorem's constant 1030 is untouched:

**REPAIR 1 — MINOR FIX (R4, eq. (14), tail): dropped factor e^{a/2}.** In (14), e^{a c_0/2} with c_0 = 1+(i+1)/√2 equals e^{a/2}·e^{a/(2√2)}·(e^{a/(2√2)})^i, but (14) keeps only e^{a/(2√2)}e^{ex'} — the factor e^{a/2} = 1.74224 is missing from the first summand, so the written derivation does not prove (14). The stated (14) is nevertheless still *true*: it exceeds the exact tail box-sum by ~10^62 at k=99 (exact 5.0074e-79 vs stated bound 9.3901e-18; step_04) because the CX6 (16) sum-trick slack dwarfs e^{a/2}. Fix: replace e^{a/(2√2)} by e^{a(1+1/√2)/2} in (14)'s first summand (and in `step_15.py` tail_sum1). Consequence: C_raw → 1025.88985257 (step_07), margin 4.1101 (was 4.1102). Theorem unchanged.

**REPAIR 2 — SUPPLEMENT (R8, mid/mid2 monotonicity).** §6 gives complete pairing arguments for q2/tail/T2/T3J, but for mid and mid2 it only checks the comparison at k=99, asserts "the ratio only improves for larger k", and certifies on a grid k∈[99,1200). A grid alone is a GAP for an all-k claim. I closed it (step_06): each comparison reduces to a linear-in-k inequality with arb-certified coefficients — e.g. mid: (3π/4)(k+1) > (√(2π)/2+0.05)²(k−2) since 2.3561945 > 1.6986277; mid2: 10.995574 > 1.9707222 — valid for ALL k ≥ 99, corroborated by a certified arb grid on [99,3000] (0 violations) and by totalC(k) decreasing 1025.89@99 → 522.43@10^6 (step_09). With this supplement R8 is analytic for all k ≥ 99; the grid is redundant.

**Optional hygiene:** (i) `step_15.py` mid2 uses the finite sum Σ₁₃⁴⁹ while proof (13) has Σ_{i≥13} (the version R8 monotonizes); infinite sum is larger by relative ~1e-49, no printed digit changes — align for literal consistency. (ii) `step_18b.py`'s h() uses a(k/2−12) vs proof's a(k/2−13) (stronger test; superseded by my analytic proof h(k)>0 ∀k≥93). Cosmetic. (iii) `rst1_common.py` is a superseded prototype (requires r_* ≥ 2c_0, different mid formula; not used by the certified chain). Repo note only.

## Lemma-by-lemma verdict table

| # | Lemma | Verdict | Key evidence (all printed) |
|---|---|---|---|
| L0 | CX6 inputs (decomposition (3)=CX6(8); \|A_k\| (4)=CX6(12); \|H\|≤2, √(π/2)x^{−3/2}; \|K\|≤2e^{−ar}, \|K−E\|≤2e^{−3ar}; CX6 (17),(18),(20); \|J_k\|≤2^{−(N−1)}) | **CORRECT** (used exactly as proved) | statements match verbatim; CX6R verdict PROVED covers (12),(13),(18),(20),\|J_k\|; step_17 tail = CX6(20) at R=200 |
| R1 | modulus-exponent bound | **CORRECT** | step_02: Padé g′(u)=u²/((1+u)(2+u)²) (sympy); \|s+q\|²−(ρ+q/√2)² identity exact, min grid 0.79289>0; (5) margin min 1.04e-13>0 for j≥1 (j=0 equality) |
| R2 | Gaussian saddle lemma (r_λ>c_0) | **CORRECT** | step_03: all saddle identities exact (r_m, Ẽ(r_m)=−Σ+λc_0/2, 2λr_λ=Σ, Ẽ(r_b)=−5Σ/4+λc_0/2, Ẽ′(r_b)=3λ, c_2, concavity); hypothesis at every use: q2 min r_a/c_0=4.63, mid 2.63 (i=12), T2 2.68, all ↑ in k; exact arb quadrature vs (11): worst ratio at uses 0.00165, worst anywhere in scan 0.0034, always <1 |
| R3 | universal box lemma | **CORRECT** | step_03: ≤1.4e-4 at uses; global worst 0.075 at (j,q)=(1,1) (outside uses), always <1 |
| R4 | T_1 assembly + split | **MINOR FIX** (e^{a/2} in (14)) | step_04: split covers i=0..k exactly once; mid2 (13) valid (exact/bound 0.438@99→0.143@1000, decreasing); tail (14) slip (see Repair 1); step_12: full exact-F assembly at k=99 (101 independent quadratures): exact 2.18033152801e-12 < bound 1.31955937771e-9 (ratio 0.00165), every group exact<bound |
| R5 | convexity Q_i | **CORRECT** (fully analytic) | step_05: ℓ″(i)=√π[u^{−1/2}+(k−2i−1)²/(4u^{3/2})]>0 (sympy diff=0); h(99)=12.78289481; h(k) ≥ c_1k−√(14πk), c_1=√π/2−π/16=0.689877384603 ⇒ h(k)>0 ∀ real k≥93 (k_crit=92.4132) — grid redundant; exhaustive Q_i≤Q_13 on [99,3000]: 0 violations |
| R6 | T_2 bound | **CORRECT** | re-derived from \|K−E\|≤2e^{−3ar}; r_{3a}=4.5834978>c_0=1.7071068 at k=99, increasing; Σ=√(3πk) subdominant |
| R7 | T_3, J_k | **CORRECT** | CX6 (18) and \|J_k\|≤2^{−(N−1)} verbatim; T3J piece 1.213888250525348e-18 matches table |
| R8 | monotonicity k≥99 | **CORRECT with SUPPLEMENT** (Repair 2) | step_06: certified arb grid [99,3000] max P(k+1)/P(k): q2 0.99998, mid 0.99484, mid2 0.96339, tail 0.69083, T2 0.99484, T3J 0.56941; analytic closures for every piece (pairing + certified linear-in-k inequalities) ∀k≥99 |
| R9 | constant C_raw | **CORRECT** | step_07 (independent arb, 90 dps): q2=1012.64180233306093617442, mid=2.61805771696862494774434, T2=10.6299604664904913310377, C_raw=1025.8898389139783896816312±1.6e-23 — digit-for-digit match with proof; rigorous variant (Repair 1 + infinite mid2 sum) C=1025.8898525699464901074637; margin to 1030 ≥ 4.1101 both |
| R10 | sanity k∈[99,400] | **CORRECT** | step_08: verify_rst1.py rerun PASS (240 s): all k pass, max \|rest\|/target = 0.21501950±1.7e-9 at k=106 (matches); step_12: TSV complete on [0,400] |
| R11 | item 1 (numerical) | **NUMERICAL — corroborated** (honestly labeled) | step_11: A_k^[2] + saddle-equation identities PROVED (sympy diff=0; 50-pt mpmath max rel diff 2e-81); rest/R_k = 0.9997061, 1.002624, 1.000035, 0.9998113, 0.9999891, 0.9999727, 0.9999883 (k=100…400, matches); \|g_k\|Ne^{θ₂} = 1.9895026, 1.8650556, 1.8130177, 1.7904286 → √π=1.7724539; arg−θ₂ = 1.5109564…1.5654485 → π/2; ReΦ residual → 0; subdominance q3: 9.83e-4/4.76e-5/4.05e-6, K−E: 4.50e-3/2.49e-4/3.37e-5 (match) |
| Cor | C_3 = 1,000,668 | **CORRECT** | step_10: δ=√(2π)−√π=0.73417442372548447512; g(N)=1030N^{3/2}e^{−δ√N} decreasing for N>16.697 (3/δ=4.08622243305); g(100)=667.35699188314026003 < 667.357 (certified upper); integer grid (100,2000] max ratio 0.97854<1; C_2+g(100)=1000667.3569918831403 < 1000668 |

## Worst-case analysis (task item 4)

R2 exact/bound at proof uses is **worst at k=99** (q2 0.00165, T2 3.9e-4) and *decreases* with k: 0.00101@400, 6.8e-4@10³, 2.3e-4@10⁴, 7.6e-5@10⁵ (step_09, mpmath exact quadrature cross-checked against arb-certified exacts to ≤6e-7). R3: 1.4e-4 at uses. Global worst over a wide (j,q,λ) scan: R2 0.0034 at (15,2,a/2), R3 0.075 at (1,1,a) — both far from the proof's uses; neither lemma ever fails. The assembled normalized bound totalC(k) is dominated by k=99 and decreases strictly to 522.43 at k=10⁶.

## Digit counts / residuals (certification quality)

- C_raw: arb ball radius ±1.6e-23 on 1025.8898389139783896816312 (≈22 certified decimals); margin to 1030 = 4.1101.
- Sanity: max ratio 0.21501950±1.7e-9 (9 digits); R2/R3 ratios ±1e-11 (arb, 110 dps).
- Independent assembly exact/bound residual at k=99: 0.0016523179 (T_1), 0.00039091 (T_2).

## Files that produced this (work dir `new_results/RST1R_referee_rest/`)

- `step_01.py` recon; `step_02.py` R1 (Padé sympy + arb stress); `step_03.py` R2/R3 algebra, hypotheses, exact-vs-bound scan; `step_04.py` R4 audit (coverage, mid2/tail exact sums, e^{a/2} slip); `step_05.py` R5 analytic closure; `step_06.py` R8 certified grid [99,3000] + all-k analytic closures; `step_07.py` R9 independent constant; `step_08` = verify_rst1.py full rerun → `verify_rst1_rerun.log` (PASS, 240 s); `step_09.py` large-k stress (k=10³,10⁴,10⁵; totalC to 10⁶); `step_10.py` corollary; `step_11.py` R11; `step_12.py` full exact-F assembly at k=99 + TSV completeness. `progress.log`: one line per finished step. `result.md`: this report (kept current throughout).

**Bottom line:** with Repair 1 (e^{a/2} in (14)) and Repair 2 (the mid/mid2 all-k monotonicity closures I supplied) applied, the proof is complete and the theorem stands exactly as stated: |M_k − M_k^main| ≤ 1030·N^{−1/2}e^{−√(2πN)} for all N ≥ 100, with C_raw = 1025.8898389… certified; the corollary C_3 = 1,000,668 is certified; item 1 is numerical as labeled and fully corroborated.