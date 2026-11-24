# SRF3R — Adversarial referee verdict on SRF3 (octahedron end)

**Subject:** `new_results/SRF3_octahedron_radius_proof/` (proof.md, result.md, verify_srf3.py, `SRF3_hi_coeffs_K120.json`), against the independent `SRF1_radii_3d/SRF1_n3_K60.json` and the 2D template (`paper/ieee/main.tex`, thm:zero). **Every number below was printed by my own code** (`verify_srf3r.py`, 136 s, all PASS; plus `step_21–61.py`). Notation: z₀ = (−log2 + iπ)/2, R = |z₀| = 1.60857525585590493, θ = arg z₀ = 1.7879528081463, bₖ = aₖ/k!.

## 1. Verdicts

| Claim | Verdict | Reason |
|---|---|---|
| **A** — S(p) = 8(1−τ)²T(τ) + square/Laplace/logit forms | **PROVED** | Algebra + all three Jacobians exact; my direct quadrature (first-octant **graph** and **radial-over-sphere**, not via T) agrees with 8(1−τ)²T to **50.5–67.5 digits** at d = ±0.01, ±0.03; forms agree 24–45+ digits |
| **B** — nearest Borel singularities exactly z₀, z̄₀, square-root | **PROVED AFTER REPAIRS** | §3.2's contour claim is false for moderate x (repair: V-split, written below); zeros, corner min, Hankel constant, square-root nature all independently verified; two uniformity steps remain proof-sketch+ (admitted by SRF3, same status as the paper's thm:zero †) |
| **C** — radius zero, lim sup \|bₖ\|^1/k = 2/√(π²+log²2) = 0.6216681478591508 | **PROVED AFTER REPAIRS** | §3.2 GAP repaired (V-split); §3.3 R2-split defective (H₀ off by eʷ; "Ψ-part analytic in \|z\|<2R" false as stated) — non-cancellation re-proved directly; §3.4 correct, trig lemma verified to k = 2·10⁵ |
| **D** — growth law (τ and d series) | **PROVED AFTER REPAIRS** | τ-law and the 2nd d-display correct, confirmed to **9–10 digits** by my own fits; **1st d-display sign error CONFIRMED** (planner's −1.54 reproduced exactly); β = −1.500 sharp; Cτ err 1.4e−9; C_d err 1.8e−10 |
| **E** — Borel summability along ℝ₊ (NUMERICAL) | Confirmed, strengthened | My Padé–Borel [55/55] vs my own quadrature: **46.0, 47.3, 46.5, 34.7 digits** at d = 0.01, −0.01, 0.03, 0.1 (SRF3: 37.0, 15.2, 6.3) |

## 2. Sign convention (Referee task 2) — RESOLVED

- **SRF1's `tau`/`d` are bit-identical to SRF3's `eta`/`d` for k ≤ 60** (step_21): same variable, same series.
- Raw ratios (printed): bₖ(d)/(k^{−3/2} Im[z₀^{1/2−k}]) for k = 52…60 = **−1.525, −1.599, −1.537, −1.128, −1.547, −1.500, −1.559, −1.525, −1.585**, oscillating about **−√(3π/4) = −1.53499** — the planner's −1.54 reproduced. With Im[z̄₀^{1/2−k}]: +1.52…+1.60. τ-ratios oscillate about −√(3π/8) = −1.08540 as claimed.
- **Resolution:** Claim D's first d-display is wrong by one sign; the mixing factor e^{−z̄₀} = √2 i (printed: 1.414213562 i) is correct; proof.md's *second* d-display is correct. The correct first display:
  bₖ(d) ~ **−√(3π/4) k^{−3/2} Im[z₀^{1/2−k}] = +√(3π/4) k^{−3/2} Im[z̄₀^{1/2−k}] = +√(3πR/4) R^{−k}k^{−3/2} sin((k−½)θ)**, confirmed by my fit C₂ = 1.94682206573 vs theory 1.9468220655487 (9.7 digits), C₁ → 0.
- Derivation check (mine): bₖ(d)/k!-mixing gives bₖ(d) ~ √(3πR/8)·Re[e^{i(kθ+φ)}e^{−z̄₀}]·R^{−k}k^{−3/2}, φ = π−θ/2, e^{−z̄₀} = √2 i ⟹ phase shift **+π/2** (cos(kθ+·) convention) ⟹ the sin((k−½)θ) law with +√(3πR/4). SRF3's "φ_d−φ_τ = 3π/2" is the same law in the atan2(B,A)/cos(kθ−·) convention; their "φ_τ = π−θ/2 = −2.2476163" is a sign typo (π−θ/2 = **+**2.2476162; −2.2476 is −(π−θ/2)).

## 3. Line-by-line marks (proof.md)

| Step | Mark | Detail / repair |
|---|---|---|
| §1 Thm A + 3 forms | **CORRECT** | Verified algebraically and to 50–67 digits by two independent parametrizations |
| §2 Borel contour (1/2πi)∮f e^{z/τ}dτ/τ | **CORRECT** | Residue check redone |
| §2 Hankel B_f ~ −iE(τ*)²/(2√π)z^{−3/2}e^{z/τ*} | **CORRECT** (up to stated determination) | Re-derived exactly via ₁F₁(−1/2;1;w) → −(−τ*)^{1/2}(τ*)^{3/2}/(2√π) z^{−3/2}e^{z/τ*} |
| §2/§4.1 zeros τ* = (z₀/x)(1 **+** log(s(1−s))/(2x)+…) | **MINOR FIX** | Correct sign is **−**: τ* = (z₀/x)(1 − log(s(1−s))/(2x)+O(x⁻²)), so \|τ*\| > R/x, min at s = ½ (verified: ratios 1.163, 1.131, 1.074, 1.036 at x = 5, 6, 10, 20 ≈ 1+log2/x). The 1/τ* display and everything downstream use the correct sign |
| §3.1 interior piece | **MINOR FIX** | Integrand is *not* entire in τ for t_i ≥ δ (Σtᵢ^{2τ} has zeros, e.g. 2t₁^{2τ} = −t₃^{2τ} at \|τ\| = R/log(t₁/t₃)); but zeros are uniformly at \|τ\| ≥ R/log((1−2δ)/δ) > 0 (same fiber equation), so coefficients grow only geometrically ⟹ entire Borel contribution. Conclusion stands |
| §3.1 corner analysis; "\|z\| = R only on edges" | **CORRECT (numerically decisive); uniformity proof-sketch+** | My all-branch root search of e^{−2zw}+e^{−2z(1−w)} = −1, 13 w-values: min \|z\|/R = 2.0 (w=½), 1.827 (0.3), 1.178 (0.1), 1.081 (0.05), 1.0076 (0.005) — **≥ 1 everywhere, → 1 only at w → 0,1 (the edges), = 2 at the fiber midpoint**; actual-Q corner scan: min \|τ*\|U/R = 0.81 (U=3), 0.93 (U=5), 1.006 (U=8), 1.044 (U=12). Double roots only at w = 1 (solving F = F_z = 0 forces w = 1), so on compact w-sets roots are simple and continuous — Rouché closes the repair |
| §3.1 kernel \|k_j\| ≤ (C(V)e^{−V})^j/j!, C(V) = 2e^V\|log(1−e^{−V})\| | **CORRECT (exact)** | k_j = L^j/j!; C(V) decreasing 4.1718 → 2 (printed at V = 0.2, 1, 5, 20) |
| §3.2 contour "inside all zeros of Q" | **GAP (false as written); repair written** | min_{s∈[δ,1−δ]} \|τ*\|x/R = **0.85 (x=2), 0.72 (x≈1.2–1.7)**, 0.61 at s=0.015 — zeros are INSIDE \|τ\| = R/x for moderate x, *within the strip*. SRF3's step_06 sampled only s ∈ {0.5,0.25,0.1} and hit a wrong Newton branch at u = 0.1 (printed 0.0056), hiding the failure. **Repair:** split at V₀(ε): for x ≥ V₀ all zeros lie outside \|τ\| = (R−ε)/x (min ratio ≥ 1.03 at x ≥ 20, attained at s = ½); the bounded-V part contributes an **entire** function to B_T (uniform zero-free disc by compactness, Q(0) = 3 ≠ 0). Corners: same split in U ≥ U₀. Limsup ≤ 1/R intact |
| §3.2 "\|Q\| ≤ 3u^{−2r}", "u^{−r}4^{2r}" | **MINOR FIX** | Holds only if s(1−s) ≥ u. Correct: \|f\| ≤ √2[u(1−u)]^r + [(1−u)²s(1−s)]^r ≤ √2+1 (Re τ = −r); (δ(1−δ))^{−r} for 4^r (Re τ = +r). Max on the circle at σ = ±r valid by log-convexity. Constants only |
| §3.2 "corner/kernel 2R with room to spare" | **MINOR** | Corners attain R/U edge-ward (2R only at w = ½); harmless |
| §3.3 H₀(w) = e^w√(1+2e^{2w}) "strict V→∞ limit" | **WRONG** | Printed: Φ(w/V,V) → √(1+2e^{2w}) = e^w√(2+e^{−2w}); ratio Φ/H₀ = 0.7499 (w=0.3), 0.5005 (w=0.7) = e^{−w} |
| §3.3 "Ψ-part analytic in \|z\| < 2R" | **WRONG as stated; conclusion re-proved** | The Ψ-part keeps z₀-strength (πi/8 − 1)i√2 z₀^{1/2} ≠ 0: per-edge strength (πi/4)e^{z₀}z₀^{1/2} vs Ĥ-strength i√2 z₀^{1/2}, ratio printed = **πi/8 exactly** (0.39269908 i to 40 digits) — the s-integral's B(3/2,3/2) factor. **Repair (the correct non-cancellation):** §4.1 computes the total strength D ≠ 0 directly; bounded-V part is entire; higher Hankel terms ∝ (1−z/z₀)^{3/2,5/2,…} cannot cancel the (1−z/z₀)^{1/2}; next ladder zeros at ½√(log²2+9π²) = **4.7251 > 2R = 3.2172** (printed). Combinatorial steps (n−j)!/n! ≤ 1/j!, (j+1)^{n−j+1} ≥ 2^{n−j+1} verified n < 60 |
| §3.4 lower bound | **CORRECT** | 3-consecutive-k lemma: 0 violations to k = 2·10⁵; weakest window max\|cos\| = 0.7796 ≥ sin((π−θ)/2) = 0.62632 (good-gap θ > bad-gap π−θ; no irrationality needed). Printed \|bₖ\|^{1/k}R ≈ 0.93 at k = 120 with the expected −(3/2)(log k)/k approach. τ↔d mixing identity worst rel err **2e−69** (k ≤ 120) |
| §4.1 singular-part computation | **CORRECT** | E → 2e^{z₀}√x (Q'(τ*) = 4x+O(1), u^{2τ*} → −2, e^{τ*M} → e^{z₀} = i/√2 printed); e^{z/τ*} = e^{zx/z₀}[s(1−s)]^{z/(2z₀)}(1+o(1)); B(3/2,3/2) = π/8; Γ(−½) = −2√π; per edge (πi/4)e^{z₀}z₀^{1/2}(1−z/z₀)^{1/2}; ×3 edges ×2/√3; Borel[τg] = ∫Borel[g] so primitives raise branch order; strength (√3πi/2)e^{z₀} = −π√(3/8) in the tracked determination; sign of D fixed by fit (admitted open phase step) |
| Thm D τ-displays | **CORRECT** | C₁ fit −1.3766110829 vs −√(3πR/8) = −1.3766110843131 (9 digits) |
| Thm D first d-display | **WRONG SIGN** (see §2 above) |
| Thm D second d-display | **CORRECT** | C₂ fit err 1.8e−10 (9.7 digits) |
| §4.2 constants; §4.3 invariance | **CORRECT** | C = (3π)^{1/2}(π²+log²2)^{1/4}/4 ✓; T(0) = √3/2; my w-fit C_w/C_τ = 0.9999999978; w₂..w₆ = paper to 2e−50; d²,d³,d⁴ to 1e−71 |
| §4.4 table | Reproduced independently | β = −1.500 (resid 2.6e−14 vs 1e−8 at ±0.002); R, θ profiles sharp (×10⁵ resid at ±1e−5 rel); vertex-angle amplitude ratio **5e−13** (< their 2e−3); phase lines have the two sign-convention typos noted |
| §5 (E) | NUMERICAL, confirmed | [55/55]: 46.0/47.3/46.5/34.7 digits at d = 0.01/−0.01/0.03/0.1 |

## 4. The two admitted open uniformity steps — referee's repair statements

1. **Corner decomposition (§3.1).** Needed: on compact w-sets [w₀, 1−w₀], min \|z(w)\| ≥ R(1+δ″) (numerically: min at endpoints, e.g. 1.081R at w₀ = 0.05); roots simple (double root ⟹ w = 1, proved); \|Q(τ,s,u) − (fiber)\| = o(1) on \|τ\| = R(1+δ″/2)/U uniformly as U → ∞ (Rouché); glue to edge strips (δ′ < δ) by a smooth partition. Status: routine, unwritten; every input verified numerically here.
2. **Hankel transfer phases (§4).** Needed: the per-(s,u) Hankel asymptotic B_f ~ −iE(τ*)²/(2√π)z^{−3/2}e^{z/τ*} with error uniform in the tail, and analytic tracking of the (τ−τ*)^{1/2} determination (currently D's sign rests on the 1e−8 cosine fit, which I reproduced to 9 digits). Status: proof-sketch+ (paper's thm:zero † level).

## 5. verify_srf3r.py (136 s < 5 min, all PASS) — printed output

```
R1 Theorem A: graph & radial vs 8(1-tau)^2 T, worst digits (d=+-0.01): 50.5
R2a tail min|tau*|x/R (x=6,12; need>=1): [1.1311, 1.0613]
R2b moderate-x min|tau*|x/R (x=2; S3.2 contour fails => V-split repair): [0.8504]
R2c corner fiber min|z|/R (w=0.5,0.3,0.1,0.05,0.02): [2.0, 1.8271, 1.1779, 1.0814, 1.031]
R3a SRF1 d,tau == SRF3 d,eta (k<=60): True
R3b d-ratio with Im[z0^(1/2-k)], k=56..60 (theory -1.53499): [-1.547, -1.5002, -1.5587, -1.5248, -1.5849]
R3c tau-law C1 (theory -1.3766110843131): -1.3766111099966
R3d d-law C2 (theory 1.9468220655487); C1~0: ('1.9468220671232', '4.282e-8')
R4a beta profile: sharp min at -1.5 (2.5e-10 vs 1e-6 at +-0.01)
R4b R profile resid at R*(1+[-1e-5,0,+1e-5]): ['2.7e-06', '3.3e-08', '2.8e-06']
R5a tau->d mixing identity worst rel err: 1.97e-69
R5b d^2,d^3,d^4 vs paper thm:oct: ['1.13e-71', '3.4e-72', '6.34e-71']
R5c 3-consecutive-k lemma violations (k<=20000): 0
R6 Pade-Borel[45/45] at d=0.01 vs T-based value (digits): inf
```

## 6. Evidence table (residuals / digit counts, all printed by my code)

| Check | Result | Where |
|---|---|---|
| Thm A: graph & radial vs 8(1−τ)²T, d = ±0.01, ±0.03 | 55.9–67.5 digits (dps 55); 50.5 in verifier | step_23, step_58, R1 |
| T: square/Laplace/logit vs simplex | exact at 45 dps / 24 digits | step_24 |
| zeros tail vs moderate-x | ≥1.03 (x≥5) / 0.61–0.85 (x≤3) | step_26/27, R2 |
| corner fiber (all branches) + actual-Q corner | min \|z\| ≥ R, 2R at w=½; ≥1 for U≥8 | step_29, step_31, R2c |
| C(V) monotone 4.1718→2; combinatorial inequalities | exact | step_34, step_35 |
| Φ(w/V,V) limit = √(1+2e^{2w}) (§3.3's H₀ wrong by e^w); strength ratio = πi/8 | 40 digits | step_33 |
| d-sign resolution; Cτ fit err 1.4e−9; C_d fit err 1.8e−10; C₁(d)→0 | 9–10 digits | step_21, step_39, step_41, R3 |
| β sharp at −1.5; R, θ sharp at theory; vertex amplitude 5e−13 | profiles printed | step_36–41, R4 |
| mixing identity 2e−69; d²d³d⁴ 1e−71; w₂..w₆ 2e−50; C_w/C_τ = 0.9999999978 | exact/low-order | step_42, step_48, R5 |
| limsup: \|bₖ\|^{1/k}R ≈ 0.93→1; 3-consec-k lemma 0 violations to 2·10⁵ | printed | step_44, R5c |
| Padé–Borel [55/55] vs own quadrature: 46.0/47.3/46.5/34.7 digits | d = 0.01/−0.01/0.03/0.1 | step_43, R6 |

## 7. Files that produced this report (work dir `new_results/SRF3R_referee_octahedron/`)

- **`verify_srf3r.py`** — consolidated verifier (source in step_59.py; run output above, 136 s).
- Key steps: step_21 (sign ratios, file identity), step_23 (direct quadrature 55–67 digits), step_24 (equivalent forms), step_26/27 (zero scans: GAP evidence + tail), step_29 (corner fiber all-branch), step_31 (actual-Q corner), step_33 (Φ limit + πi/8 strength ratio), step_34/35 (kernel bound, combinatorics), step_36/38/39/40/41 (fits: C, β, R, θ, vertex amplitude), step_42 (exact low order, w-series), step_43 (Padé–Borel), step_44 (3-consecutive-k lemma + limsup), step_48 (mixing identity, phase fits), step_49/50–58 (phase-convention resolution; float-vs-decimal p pitfall identified — the apparent 17-digit graph/T mismatch was float `1.01` vs decimal `mpf('0.01')`, explained exactly by dS/dp ≈ 6.9), step_55–57 (independent tanh-sinh engine validation).
- `result.md` (full referee report) and `progress.log` (one line per step) updated.

**Bottom line:** A is PROVED. B, C, D are PROVED AFTER REPAIRS: (i) Claim D's first d-display sign flip (planner's finding confirmed: −√(3π/4), not +); (ii) §3.2's contour claim replaced by the V-split (bounded-V ⟹ entire; tail zeros verified outside (R−ε)/x); (iii) §3.3's split replaced by the direct non-cancellation argument (H₀'s e^w slip and the πi/8 strength mismatch documented); (iv) sign typos in the τ* display and the two phase-table lines. Every constant in the growth law (A, β = −3/2, Cτ, C_d = √2 Cτ, both phases) is independently confirmed to 9–10 digits, and the two admitted uniformity gaps remain at the paper's thm:zero-† level with concrete repair paths and decisive numerical support.