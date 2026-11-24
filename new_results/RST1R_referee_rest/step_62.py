result = '''# RST1R referee report — sharp rest exponent sqrt(2 pi) in the squircle decay law (Kimi RST1)

Under review: `research/08_theory/RST1_rest_exponent/{proof.md, result.md, verify_rst1.py,
rst1_common.py, rst1_common2.py, step_15.py, step_17.py, step_18b.py}`.
Claim: for every N = k+1 >= 100, |M_k - M_k^main| <= 1030 N^{-1/2} e^{-sqrt(2 pi N)}
(C_raw = 1025.8898389..., arb), plus corollary C_3 = 1,000,668, plus item 1 (numerical).
Referee work dir: `research/08_theory/RST1R_referee_rest/`. Every number below was printed
by the step_NN.py files in this directory (logs on stdout; verify_rst1 rerun in verify_rst1_rerun.log).

## OVERALL VERDICT: **PROVED AFTER REPAIRS**

The theorem is correct and every load-bearing lemma re-derives independently. Two repairs are
required (both far below the certified margin; the constant 1030 is untouched):

1. **MINOR FIX (R4, eq. (14) tail): dropped factor e^{a/2}.** In (14) the box prefactor is
   e^{a c_0/2} with c_0 = 1 + (i+1)/sqrt2, i.e. e^{a/2} * e^{a/(2 sqrt2)} * (e^{a/(2 sqrt2)})^i.
   The displayed (14) keeps only e^{a/(2 sqrt2)} e^{e x'} — the factor e^{a/2} = 1.74224 is
   missing from the first summand. The derivation as written therefore does not prove (14).
   Numerically the stated (14) still exceeds the exact tail box-sum by a factor ~ 5e-62 / 1
   at k = 99 (step_04: exact 5.0074e-79 vs stated 9.3901e-18) because the CX6 (16) sum-trick
   slack dwarfs e^{a/2}; so the final inequality is still true. Repair: replace
   e^{a/(2 sqrt2)} by e^{a(1+1/sqrt2)/2} in (14)'s first summand (and in step_15.py tail_sum1).
   Effect on C_raw: +1.3656e-5 -> C_raw = 1025.88985257 (step_07), margin to 1030 = 4.1101
   (was 4.1102). Theorem unchanged.
2. **SUPPLEMENT (R8, mid/mid2 monotonicity):** proof.md SS6 proves q2/tail/T2/T3J termwise by
   complete pairing arguments, but for mid and mid2 it only checks the comparison at k = 99 and
   asserts "the ratio only improves for larger k", with grid certification on [99,1200). A grid
   alone is a GAP for an all-k claim. Repair (step_06, supplied by referee): the comparisons
   reduce to linear-in-k inequalities with arb-certified coefficients, e.g.
   (3 pi/4)(k+1) > (sqrt(2 pi)/2 + 0.05)^2 (k-2) for mid, using 1/(2N) <= 0.05/sqrt N;
   these hold for ALL k >= 99 (and mostly k >= 1). Plus a certified arb grid on [99,3000]
   (0 violations). R8 then holds analytically for all k >= 99; the grid is redundant.

Optional hygiene: (i) step_15.py computes mid2 with the finite sum xsum = sum_{13}^{49}, while
proof eq. (13) has sum_{i>=13} (the version R8 monotonizes); the infinite sum is larger by a
relative ~1e-49 — no printed digit of C_raw changes, but literal consistency is preferable.
(ii) step_18b.py's h() uses a(k/2-12)/(2 sqrt2) where proof.md and verify_rst1.py use
a(k/2-13)/(2 sqrt2); the step_18b version is a *stronger* test, superseded by the analytic
proof of h(k) > 0 for all k >= 93 (step_05). Cosmetic.
(iii) rst1_common.py is a superseded prototype (its FH_gauss requires r_* >= 2 c0 and uses a
relaxed mid formula; its T1_bound boxes all i>=2). Not used by verify_rst1.py/step_15/step_18b;
no impact on the theorem. Repo hygiene note only.

## Lemma-by-lemma verdict table

| # | Lemma | Verdict | Evidence (all printed) |
|---|---|---|---|
| L0 | CX6 inputs (decomposition (3)=CX6(8); |A_k| (4)=CX6(12); |H|<=2, |H|<=sqrt(pi/2)x^{-3/2}; |K|<=2e^{-ar}; |K-E|<=2e^{-3ar} (13); CX6 (17),(18),(20); |J_k|<=2^{-(N-1)}) | **CORRECT** (used exactly as proved) | read both proofs; CX6R verdict PROVED covers (12),(13),(18),(20),|J_k|; step_17 tail = CX6(20) at R=200 verbatim |
| R1 | modulus-exponent bound | **CORRECT** | step_02: Pade g'(u)=u^2/((1+u)(2+u)^2) (sympy); |s+q|^2-(rho+q/sqrt2)^2 = q+q^2/2+sqrt2 q(r-rho) verified exactly, min grid value 0.79289>0; (5) margin min 1.04e-13>0 for j>=1 (j=0 equality); rho<=r+1/2 OK |
| R2 | Gaussian saddle lemma (r_lambda>c_0) | **CORRECT** | step_03: all saddle identities (r_m, E~(r_m)=-S+lam c0/2, 2 lam r_lam=S, E~(r_b)=-5S/4+lam c0/2, E~'(r_b)=3 lam, c2, concavity) exact; hypothesis at every use: q2 min r_a/c0=4.63, mid min 2.63 (i=12,k=99), T2 2.68, all increasing in k; exact arb quadrature vs (11): worst ratio at proof uses 0.00165, worst over extra scan 0.0034, always <1 |
| R3 | universal box lemma | **CORRECT** | step_03: exact/bound <= 1.4e-4 at proof uses; worst over wide scan 0.075 at (j,q)=(1,1) (outside uses), always <1 |
| R4 | T_1 assembly + split | **MINOR FIX** (e^{a/2} in (14)) | step_04: split covers i=0..k exactly once; mid2 (13) valid (exact/bound 0.438@99, 0.224@400, 0.143@1000, decreasing); tail (14) missing e^{a/2} yet still > exact sum by ~1e62; per-term prefactor bounds OK incl. k=27 corner; step_12 full exact-F assembly at k=99 (101 independent quadratures): exact 2.18033152801e-12 < bound 1.31955937771e-9 (ratio 0.00165) |
| R5 | convexity Q_i | **CORRECT** (closed analytically) | step_05: ell''(i)=sqrt(pi)[u^{-1/2}+(k-2i-1)^2/(4u^{3/2})]>0 (sympy diff=0); h(99)=12.78289481; h(k) >= c1 k - sqrt(14 pi k), c1=sqrt(pi)/2-pi/16=0.689877384603, so h(k)>0 for ALL real k>=93 (kcrit=92.4132) — grid checks redundant; exhaustive Q_i<=Q_13 k in [99,3000]: 0 violations; ledger's Q_i<=Q_2 (unused in assembly) holds for k>=21 |
| R6 | T_2 bound | **CORRECT** | re-derived (2/pi)|T_2| <= (4/pi) F_{H,3a}(k,1) from |K-E|<=2e^{-3ar}; r_{3a}=4.5834978>c_0=1.7071068 at k=99, increasing; S=sqrt(3 pi k) subdominant; step_03/step_09 ratios |
| R7 | T_3, J_k | **CORRECT** | CX6 (18) (2/(pi k))3^{-k} and |J_k|<=2^{-(N-1)} used verbatim; step_07 T3J piece = 1.213888250525348e-18 matches proof table |
| R8 | monotonicity k>=99 | **CORRECT with SUPPLEMENT** (mid/mid2 closure supplied) | step_06: certified arb grid [99,3000], max P(k+1)/P(k): q2 0.99998, mid 0.99484, mid2 0.96339, tail 0.69083, T2 0.99484, T3J 0.56941; analytic closures for ALL k>=99 for every piece (pairing + certified linear-in-k inequalities, e.g. mid: 3pi/4=2.3561945 > (B+0.05)^2=1.6986277, margin at k=99 = 70.85>0; mid2: 10.995574 > 1.9707222; tail: RHS(99)=0.1403819<1/2; T3J: -0.5628<0). step_09 Part 2: totalC(k)=1025.89@99 -> 522.43@1e6 strictly decreasing |
| R9 | constant C_raw | **CORRECT** | step_07 (independent arb, 90 dps): q2=1012.64180233306093617442, mid=2.61805771696862494774434, T2=10.6299604664904913310377, C_raw=1025.8898389139783896816312 +/- 1.6e-23 — matches proof's table digit for digit; proof-rigorous variant (e^{a/2} tail + infinite mid2 sum) C=1025.8898525699464901074637; margin to 1030 >= 4.1101 both variants |
| R10 | sanity k in [99,400] | **CORRECT** | step_08: verify_rst1.py rerun PASS (240 s): all k pass, max |rest|/target = 0.21501950 +/- 1.7e-9 at k=106 (matches proof); step_12: TSV complete [0,400], no missing k |
| R11 | item 1 (numerical) | **NUMERICAL — corroborated** (labelled numerical, not proved) | step_11: A_k^[2] identity and saddle-equation polynomial PROVED (sympy diff=0, plus 50-point mpmath spot check max rel diff 2e-81); rest/R_k = 0.9997061, 1.002624, 1.000035, 0.9998113, 0.9999891, 0.9999727, 0.9999883 at k=100..400 (matches proof); |g_k|Ne^{th2} = 1.9895026, 1.8650556, 1.8130177, 1.7904286 at k=400,2000,1e4,5e4 -> sqrt pi=1.7724539; arg-th2 = 1.5109564, 1.544049, 1.5588371, 1.5654485 -> pi/2; Re-Phi residual -> 0 (-0.0425,-0.0187,-0.0083,-0.0037); subdominance q3: 9.83e-4/4.76e-5/4.05e-6, K-E: 4.50e-3/2.49e-4/3.37e-5 (match) |
| Cor | C_3 = 1,000,668 | **CORRECT** | step_10: delta=sqrt(2pi)-sqrt(pi)=0.73417442372548447512; g(N)=1030 N^{3/2}e^{-delta sqrt N} decreasing for N>16.697 (3/delta=4.08622243305); g(100)=667.35699188314026003 < 667.357 (certified upper); integer grid (100,2000] max ratio 0.97854<1; C_2+g(100)=1000667.3569918831403 < 1000668 |

## Worst-case analysis of the uniform bounds (task item 4)

- R2 (Gaussian): exact/bound at proof uses is worst at k=99 (0.00165 q2, 3.9e-4 T2) and
  *decreases* with k: 0.00101@400, 6.8e-4@1e3, 2.3e-4@1e4, 7.6e-5@1e5 (q2). Over a wider
  (j,q,lam) scan the worst ratio anywhere is 0.0034 at (15,2,a/2); the lemma never fails.
- R3 (box): worst at uses 1.4e-4; global worst in scan 0.075 at (j,q)=(1,1) — far from uses.
- Conclusion: the constant C_raw is dominated by k=99 (the norm of the bound relative to the
  true rest improves monotonically afterwards; totalC(k) decreases to 522 at k=1e6).

## Repairs required (restated)

1. (14): insert e^{a/2} in the first summand: e^{a/(2 sqrt2)} -> e^{a(1+1/sqrt2)/2}
   (and step_15.py tail_sum1 accordingly). C_raw -> 1025.88985257; margin intact.
2. R8 mid/mid2: adopt the referee's linear-in-k closures (step_06 Part B) in place of the
   k=99 comparison + grid.
3. (optional) step_15.py mid2 xsum -> infinite sum sum_{i>=13} to match (13) literally.
4. (optional) align step_18b.py h() with proof (k/2-13), and mark rst1_common.py superseded.

With 1+2 applied, the proof is complete and the theorem stands as stated with C = 1030.

## Files in this referee dir

- step_01.py recon; step_02.py R1; step_03.py R2/R3 (+hypotheses); step_04.py R4 audit
  (coverage, mid2/tail exact sums, e^{a/2}); step_05.py R5 analytic closure; step_06.py R8
  certified grid + analytic closures; step_07.py R9 independent constant; step_08 =
  verify_rst1.py full rerun (verify_rst1_rerun.log, PASS 240 s); step_09.py large-k stress
  (k=1e3..1e5 exact-F vs bounds, totalC to 1e6); step_10.py corollary; step_11.py R11;
  step_12.py full exact-F assembly at k=99 + TSV completeness.
- progress.log: one line per step. result.md: this report.
'''
open('result.md','w').write(result)
open('progress.log','a').write("step_13: result.md finalized: verdicts L0-R11+corollary; overall PROVED AFTER REPAIRS (minor fix e^{a/2} in (14); R8 mid/mid2 closure supplied; 2 optional hygiene notes).\n")
print("result.md written,", len(result), "chars")
