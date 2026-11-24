report = r'''# DCY2R — referee report on DCY2 (first correction ϖ₁ of the decay law of M_k)

**Verdict: PROVED AFTER REPAIRS.** The theorem (main term, C₂ = 10⁶) and its corollary
(full M_k, C₃ = 1.854×10¹⁰) stand. One lemma (L11, the finite-range check) was **WRONG as
coded** — a sign error — and is repaired here (re-certified: max residual 1.7512 ≪ C₂).
Several smaller gaps/typos are repaired below; none affects the constants.

## 1. The discrepancy (planner finding) — RESOLVED

The planner's quantity is the theorem's: max over N=100…401 of
|M_k + (√π/N)e^{−θ}Re[e^{iθ}(1+ϖ₁/√N)]|·N²e^{θ} = **1.744** at k=157 (N=158), full M_k.

DCY2's finite-range code (`step_06_finite.py`, `verify_dcy2.py` §3) computes
`res = abs(M - (-T1))` with `T1 = -(√π/N)e^{-θ}Re[...]`, i.e. |M_main − X| with
X = +(√π/N)e^{−θ}Re[...] — the **wrong sign** (theorem needs |M_main + X|).
Since M_main ≈ −X (paper's obs:Mk), DCY2's `res` ≈ 2|X| = 2(√π/N)e^{−θ}|Re[...]|,
and `res/bound` ≈ 2√πN|Re[...]|/10⁶, which **grows like N**: its max on [100,400] sits at
the largest N (381) and equals 1409.47 — exactly DCY2's "constant ≈ 1409" and
"max res/bound 0.00141". So DCY2's number is the artifact 2√π·381·|Re[e^{iθ}(1+ϖ₁/√381)]|,
NOT a residual. Certified (arb, outward, both signs at N=381):

    |M_main - X| N^2 e^th = 1409.4657     <- what DCY2 measured (sign error)
    |M_main + X| N^2 e^th = 1.6196351     <- the theorem's quantity

Consequence: as committed, L11 certifies a trivial statement (2|X| ≤ C₂N^{-2}e^{−θ}) and the
range N ∈ [100,400] is **not** covered (region B starts at 400). The same wrong sign appears
in `verify_dcy2.py` §3 and in the numbers quoted for L14
("residual = 0.000165, 0.000751, 0.000932" at N=100,200,400 are the same artifact).

**Repair (done here, certified):** with the correct sign, all 301 certified quadratures give
    max |M_main + X| N^2 e^th = 1.751172383 at N=116
    max |M_k    + X| N^2 e^th = 1.743787030 at N=158 (full M_k, certified table balls)
so the finite range passes with C₂ = 10⁶ (margin 5.7×10⁵), matching the planner's 1.744.

## 2. Lemma verdicts

| # | Lemma | Verdict | Notes |
|---|-------|---------|-------|
| 1 | Ray representation (CX6) | CORRECT | M_main=(2/π)Im∫_ray H E Q_k ds; code matches; table cross-check |M_k−M_main|~1e-13 |
| 2 | ϖ₁ closed form | CORRECT | independently re-derived from φ(s)=iπs/2+(N−1)log s−N log(1+s): G₁,G₂,c₁,S₁ all match; diff vs paper = 0 exactly; 40-digit match |
| 3 | Factorization, modulus, single-max, tail rates | MINOR FIX | factorization exact but for J:=(2/π)Ĩ (M_main=Im J) — §1 writes Ĩ for both (Ĩ_ray≈(π/2)T₁, so "|Ĩ−T₁|≤C₂N^{-2}e^{−θ}" as printed is false; the code consistently analyzes J). P(u) u-coefficient misprinted: √2(N−2)y²−2b_R**y²** (not y³; same typo in lamL code, harmless — min selects πy/2). P(0⁺)>0, P(1)=−2y−(2√2+π/√2)y²−2y³<0, P″<0 verified → single root; tail rates rDelta (0.514), eLamD (1.117), λ_L all verified numerically + arb |
| 4 | ℓ-boxes, c₂, Reℓ₁ | CORRECT | Reℓ₁=−y(1+b_Ry/Q₀) verified exactly (symbolic); c₂≥0.9b_Ry certified (min ratio 0.9648 on [400,1e20]) |
| 5 | Watson constants | GAP (script) / constants OK | `step_03_watson.py` FAILS its own assert (Cauchy gives M/4=0.915>0.5; interval part returns 0 via arb(lo,hi) misuse). Truth: max|ρ₃|/u³=0.0117 on (0,1], 0.053 on [1,20] → C_ρ≤0.5 and ≤2.42 both true, easily re-certified. Constants actually used are independently certified here: C_G: truth 0.746–0.779 ≤ 1 (N=400,3e3,1e5); |h_j|/((3/2)_jy^j) = 0.978/0.957/0.936 at N=400 ≤ 1+2y; M₄′ spot-checked on window grids (≤0.05) |
| 6 | Moments, model=T₁+O, the two exact identities | CORRECT | both identities re-verified exactly (diff=0); D(N) pointwise 1.657→1.565; moment-sum O(y) coeff = 0 |
| 7 | Remainder E(t) | GAP (repairable) | kept/dropped monomial algebra verified symbolically (P, drop, β₂, β₃ consistent); audit E_mom≥truth ratios 0.001–0.006 at reliable N. But: TailMom uses the one-sided formula (½)(c₂m/2)^{−s}Γ(s,z) for a two-sided integral (valid bound = 2×, exactly absorbed by slack); Mtilde misses the provable correction 2e^{a²/4c}tR^{j+1}/(j+1) (j=0 margin 3.5%, j=0 unused; j≥4 margin ≥30%). Repair: E_tailmom×2 and add the correction ⇒ worst pointwise LHS 320253→325886 (+1.8%), worst cell 390055→~397k. C₂ unchanged |
| 8 | Window-tail + single-dip (W) | CORRECT | W logic sound (W(1)<0, W″>0); the conclusion (rates valid) verified directly on grids + asymptotics |
| 9 | Region B sweep | CORRECT + MINOR FIX | 279 cells reproduced exactly, worst = 390055.4767 at (3075.23,3136.74), coverage [400,1e5] exact. Endpoint gaps ~1% (harmless): TailMom(p=0,1,2) decreasing but taken at Nhi; e^{γT³} factor at Nlo while γT³ grows ~1.2%/cell (compensated ≥1.01 by Mtilde at Nhi); assert omits +C_b=29.5 (absorbed) |
| 10 | Region C | CORRECT | cell-form LHS(1e5)=74115.35 reproduced; monotone on 923-pt dense grid [1e5,1e7] and at 1e9,1e12,1e16,1e20 (17592,14711,14090,14014 pointwise form); piecewise y-exponents ≥ y^{3/2} with log-derivative conditions at κ²=28.2 verified |
| 11 | Finite range [100,400] | **WRONG (as coded)** → REPAIRED | sign error (§1 above); re-certified: max 1.7512 ≤ C₂ |
| 12 | Part (b), C_b=29.5 | CORRECT | all six C_F constants verified: C_S≤0.596✓, C_g≤1 (truth 0.78)✓, C_S2≤2.543 (truth 1.135)✓, C_E≤3.164 (truth 0.67)✓, C_MS≤0.3555 (truth 0.29)✓, cross≤2.496 (truth 2.4951)✓; C_F=8.75≤10.155; C_b=25.38≤29.5. (§5's extra factor (1+ϖ₁/√N) is spurious — correct: |model−T₁|N²e^θ≤(π^{3/2}/2)|δ_F|/y²·(π/(2N))·... = (π^{3/2}/2)C_F) |
| 13 | Rest combination | CORRECT | N₃=25711 (boundary verified: 1.00012 vs 0.99948), h(309)=1748.9769, C₃=1.8540155e10 |
| 14 | Validation at 100,200,400 | WRONG numbers (artifact of #11); conclusion right | true scaled residuals at N=100,200,400: 0.315, 1.610, 1.335 |

## 3. Audit-script caveat (dcy2_audit.py)

Its `It` (mpmath ray quadrature, dps=50) is **unreliable for N ≳ 5·10³** (silent precision
loss; at N=1e5 it reports |It−Model| = 3.7e-250, i.e. scaled ≈ 9752, and the "E_total FAIL"
at N=1e5 is a numerical artifact). Recomputed reliably (window decomposition, H₀ restored,
cross-checked at N=400/2000 against direct quadrature to 8+ digits):

    |J - model| N^2 e^th = 0.00092 (400), 0.00055 (2000), 0.00027 (1e4),
                           8.7e-5 (1e5), 2.8e-5 (1e6)   — decreasing.

The certified bounds dominate truth by 10⁸–10⁹× at large N.

## 4. Certified constants (independently recomputed, outward-rounded)

| constant | claimed | recomputed |
|---|---|---|
| max |M_main+X|N²e^θ on [100,400] | (1409 artifact) | **1.751172383** (N=116) |
| max |M_k+X|N²e^θ on [100,400] (planner: 1.744) | — | **1.743787030** (N=158) |
| region B worst cell | 390,055 | **390,055.4767** at (3075.23,3136.74) |
| region C LHS(1e5) | 74,115 | **74,115.35233** |
| C_F (part b) | 10.155 | assembly = **8.7457**; truth 0.60 |
| C_b | 29.5 | **25.375** |
| N₃ | 25,711 | **25,711** (boundary checked) |
| h(309) | 1748.98 | **1748.976914** |
| C₃ | 1.85×10¹⁰ | **1.8540155×10¹⁰** |
| ϖ₁ | paper closed form | exact diff = 0; 40-digit match |

## 5. Required repairs

R1 (required, done): fix the sign in `step_06_finite.py` and `verify_dcy2.py` §3
   (`res = abs(M - T1)`, T1 = −X). Re-certified here: max 1.7512 ≤ C₂.
R2 (write-up): state the complex estimate for J := (2/π)Ĩ with M_main = Im J
   (the code's "Ĩ"/"model"/"target" are all J-quantities); §1's Ĩ is then the bare ray integral.
R3 (typos): P(u) u-coefficient −2b_Ry³ → −2b_Ry² (proof.md and `lamL_of`/`Pofu` code);
   §5 spurious (1+ϖ₁/√N) factor; dcy2_cert comment "|h1+(3/2)y|≤C_H1R y²" (it's |h₁|≤1.5y(1+2y)).
R4 (moment bounds): E_tailmom ×2 (provably-valid two-sided form); Mtilde correction
   2e^{a²/4c}tR^{j+1}/(j+1) (or prove the empirical domination). Impact ≤ +2% (worst cell).
R5 (sweep endpoints): TailMom(p=0,1,2) at Nlo (decreasing); e^{γT³} consistent endpoint
   (or ×1.015); EdgeInt Γ-subtraction endpoint convention. Impact ~1%.
R6 (certificates): repair `step_03_watson.py` (Cauchy M/6 + correct interval `arb` cells);
   replace the audit's `It` quadrature for N ≳ 5·10³ by the window decomposition.

## 6. Proposed paper statement

**Theorem (first correction, computer-assisted).** Put N = k+1 ≥ 100, θ = √(πN), and
ϖ₁ = [(24+12π+π²)+i(24−12π+π²)]/(48√π) = 0.84121508776758439708…−0.04501183768517361657…i.
Then, with M_k^main the ±iπ/2 branch contribution,
    |M_k^main + (√π/N)e^{−θ} Re[e^{iθ}(1+ϖ₁/√N)]| ≤ 10⁶·N^{-2}e^{−θ}      (all N ≥ 100),
and, adding the proved rest bound |M_k−M_k^main| ≤ 10,600,000·e^{−2√N},
    |M_k + (√π/N)e^{−θ} Re[e^{iθ}(1+ϖ₁/√N)]| ≤ C₃·N^{-2}e^{−θ},
C₃ = 1.86×10¹⁰ for 100 ≤ N ≤ 25,711 and C₃ = 10⁶+1 beyond. On 100 ≤ N ≤ 400 the certified
constant is ≤ 1.76. The proof: 301 certified quadratures; 279 interval cells on [400,10⁵]
(worst 3.91×10⁵); a monotone analytic tail on [10⁵,∞). In particular the O(N^{-1}) term and
the value of ϖ₁ in obs:Mk are now proved; only the rest exponent √(2π) remains numerical.

## 7. Files

`verify_dcy2r.py` (132 s): §1 ϖ₁ (sympy, exact); §2 discrepancy exhibit (N=381, both signs);
§3 repaired finite range (301 arb quadratures + table); §4 region B (their LHS_cell, 279 cells);
§5 region C; §6 part (b) constants; §7 rest; §8 tail-machinery (arb); §9 repair deltas.
`progress.log`. This report: `result.md`.
'''
open('result.md','w').write(report)
open('progress.log','a').write("step: result.md written; FINAL verdict PROVED AFTER REPAIRS\n")
print("result.md written", len(report))
