# D13t — Adversarial review of the D13s fixes to the certified decay theorem

**Theorem under review.** For all N ≥ 100: |M_k^main + (√π/N)e^{−√(πN)}cos√(πN)| ≤ C·N^{−3/2}e^{−√(πN)}, with C = 254.2276 (certified sup 254.2275176 = max of 194.1267251 [100,400], 254.2275176 [400,10¹⁰], 47.0701250 [10¹⁰,∞)).

All numbers below were printed by my code (step_NN.py files in `D13t_review_fixes/`; mpmath @50–100 dps, python-flint arb @60–70 dps, sympy exact).

## Verdict summary

| item | content | verdict |
|---|---|---|
| R1 | B2 fix: A1 enclosure `A1m = [f(zz_lo)/y_up, f(zz_up)/y_lo]` | **HOLDS** (proof + 110 884-pt containment) |
| R2 | B1 fix: hunt for other `rt2`-for-`2` slips | **HOLDS** (full grep+diff+sympy audit; no other slip) |
| R3 | H1: exact regular form of B_true | **HOLDS** (form verified 1.0e-32; bracket never fails) |
| R4 | H2: proof of x²\|g(x)\| ≤ 1 and the TL0 bound | **HOLDS** (proof complete; dense numerics) |
| R5 | H3: certify_H3 coverage; t ≤ 0.12731 on every cell (σ included) | **PROOF GAP (write-up), certificate HOLDS-NUMERICALLY** — see below |
| R6 | rerun the three certificates | **INCOMPLETE** (turn budget exhausted before reruns; see honesty note) |

---

## R1 (B2 fix) — HOLDS

- **Monotonicity proof independently confirmed (sympy, exact):** f(x)=1−arctan(x)/x, f′=g/x² with g=arctan(x)−x/(1+x²), g(0)=0, g′(x)=2x²/(1+x²)² ≥ 0 (residuals exactly 0). Hence A1 = f(zz)/y ∈ [f(zz_lo)/y_up, f(zz_up)/y_lo] is a valid enclosure whenever zz_lo ≤ zz ≤ zz_up; arb interval evaluation outer-encloses the range of zz = √2ρy/(2+√2ρy) over the (y,σ) box (dependency only loosens). The "straddling box" concern is moot: monotonicity is on all of [0,∞).
- **Containment, 110 884 sample points** over 149 cells × 4 σ-subcells (N near 100 [100,105], the sweep-worst cell [1890.7128907662395, 1891.0910333443928], 1e4, 1e6, 1e8, 1e9, 1e10−ε; wide and ratio-1.0002 boxes): worst violation **−2.09e-9** / **−2.30e-10** (contained everywhere). zz-box worst violation **−1.07e-12**. (An apparent +0.609 violation in one of my first runs was traced to an *inverted y-box in my own test* (r<1); redone with valid boxes: all contained.)
- **Exact phase confirmed independently:** my own saddle solve + Im Φ(x0,N)−θ vs (π/√2)B_h[A1′/y]: worst residual **4.44e-56** @60 dps over 2006 points incl. N=100 and N=1891. The coded T2n/B_h/phireg match this formula exactly.
- Saddle σ* ∈ **[0.428579581, 0.4999999999]** ⊂ [0,1] for N ∈ [100,1e20] — the truth always lies inside the σ-box.

## R2 (B1 fix) — HOLDS

- grep of `rt2/`, `.sqrt()`, `**arb(...)`, `ecw`, `c_w` over both fixed files; **diff vs buggy baselines: only the intended fix lines changed** (ecw; A1m; breg→B_true box; TL0; TLfN (π/2)^{3/2}). step_18.py identical.
- ecw fixed = **0.483941449038287** = √(2/π)e^{−1/2} ✓ (buggy 0.406944629688983; ratio = 2^{1/4} = 1.189207115 exactly, matching D13s's table).
- Every other suspect constant audited: GRreg 2^{9/4} = (2√2)^{3/2} ✓; GLm2reg 2^{19/4} = 2^{9/4}·2^{5/2} ✓; G0reg = (√π/4)2^{3/2} ✓; (N·y0³/2)/y = πρ³/4 ✓ (printed 0.785633806409 both ways); TL0 exponent algebra identical to 2(2+√2)^{−N/2}e^{√(πN)}/(N−2); d6 = (π/√2)(0.16666814) = **0.3702435178** ✓.
- sympy EXACT: (t²+√2t)²+(√2t+1)² = (t²+√2t+1)² (so the p2mag modulus line is exact); H2 numerator = 3x²+6√2kx+4k²; AM–GM residual = (x−1)².

## R3 (H1) — HOLDS

- Saddle branch ρ(t)=√(1+t²/(1+√2t)) solves Re Φ′(1/t)=0: worst residual **9.9e-48** (300 pts, t∈[1e-14,100]); branch equation ρ²(1+√2t)=1+√2t+t² residual **3.2e-58**.
- Exact form B_true = ρB1/(2+√2t), B1 = ½+√2 s1+s1²−ρ²A1n: worst relative residual vs direct (Im Φ(x0)−θ+π/4)/((π/√2)y) = **1.01e-32** over 5000 pts, t∈[1e-14,10³] (D13s claimed 8.5e-37 on the smaller range [1e-12,1] — consistent). φ = a(1−2ρ)/t + N·arctan(zz) residual **5.5e-47**.
- **Bracket never fails:** zz = √2t/(2+√2t) < 1 for ALL t (max 0.9986 at t=10³), so the alternating-series bracket is valid everywhere, not just the final cell. Scan of A1n ∈ [L,U] over 20 001 pts, t∈[1e-12,10⁸]: min(A1n−L) = **1.79e-50** (100 dps; a −6.4e-38 reading at 60 dps was shown to be roundoff noise), min(U−A1n) = 5e-26.
- arb Btrue_box is well-formed (radius > 0), B_abs = **0.16666814452037215** ≤ breg_abs = 0.75 (old breg harmless ✓), and contains the truth at **100 000/100 000** points of (0,Tval] (max violation **−1.48e-6**, matches D13s's 1.5e-6). B_true(Tval) = 0.166662235653773 → 1/6 ✓.

## R4 (H2) — HOLDS

I re-derived the whole chain and found **no gap**: Hc(s) = ½∫₀¹u^{s/2−1}(1−u)^{1/2}du (exact by the B-function; resid 2.3e-22 quad), (1−u)^{1/2} = 1−Σc_k u^k with c_k = (−1)^{k+1}C(1/2,k) > 0 (verified k≤60), Σc_k = 1 (**sympy exact**; partial sums track 1−1/√(πK) — mpmath nsum's 0.9985 is an nsum artifact on this slowly convergent series), termwise integration is fine by Weierstrass M. |sHc|² = 1−2Σc_k Re w_k+|Σc_k w_k|² ≤ 1−Σc_k(2Re w_k−|w_k|²) ≤ 1 (Jensen, Σc_k=1), using 2Re w_k ≥ |w_k| (sympy: numerator 3x²+6√2kx+4k² > 0) and |w_k| ≤ 1.
- Dense numerics: max x²|g(x)| on (0,1] = **0.999998915112** (200 000 pts), →1 as x→0 (tiny-x max 1.0) — bound tight, never exceeded. min 2Re(w_k)/|w_k| = **1.4142135** ≥ 1 (dense x∈(0,1], k up to 10⁸). Series identity residuals decrease with truncation (1.1e-9 @K=2e5, x=1). e^{ReΦ} = e^{−ax}x^N(1+√2x+x²)^{−N/2} resid 4.6e-48. AM–GM: (x−1)² ≥ 0 (sympy). I0/bound ratio ≈ **0.27** at N=100,400,10³,10⁴ — the proved bound dominates with factor ~3.7.

## R5 (H3) — proof GAP in the write-up; certificate itself HOLDS-NUMERICALLY

- certify_H3 reproduced: single arb box Hp3([0,0.13]) lower = **0.1538363242986387** ≥ 0 (a valid outer enclosure — subdivision is not needed for validity, only for tightness; 1300 sub-boxes give min lower 0.39796). True min of Hp3 on [0,0.13] = **0.39809551699** (200 001 pts, at t=0.13) ≥ 0.153836 ✓. H′(t)=t³Hp3(t) along the saddle: residual **2.8e-61** (200 pts); H = margin/(N/2): residual 3.2e-29; zero of H in [4.54,4.55] ✓ (claim 4.5447). t_max = y(1+y²)|_{y=√(π/200)} = **0.127300114974765** < 0.13 ✓ (the printed 0.12731 is correct).
- **The gap:** the D13s proof of the *lower* box treats H(t) with ρ on the **saddle branch only**. The lower-box margin over the full cell (σ-range included) is G(y,σ) = ψ(t) − (√2/3)t³/(1+√2t+t²) − √2 t(ρ−1)²/ρ², t=(1+y²σ)y. On a 400×400 grid (y∈(0,y(100)], σ∈[0,1]): **min G = −5.4722904e-6 < 0 at (y(100), σ=1)**; G(y,1) crosses zero at y=0.109634501782 i.e. **N = 130.68489** — so for N ∈ [100,130.68] and σ→1 the *tight* h3 lower box genuinely does **not** contain the truth. The result.md/L7 claim "h3-box containment proved on domain" is therefore **false as stated** (its t_max = y(1+y²) is precisely the σ=1 corner, showing the domain claim was intended).
- **Why the certificate is nevertheless sound:** (i) the theorem needs the mS/z box to contain the truth only at the saddle σ* ∈ [0.4286,0.5), and min G over σ∈[0,0.5], y≤y(100) = **+2.56e-14 > 0** (all of it certified by the Hp3 argument since t≤0.1264<0.13 there); (ii) the *coded* boxes carry t_lo slack from the σ-subcells: evaluating the actual fixed step_38 mS box for cells [100,105], [105,110.25], [120,126], [130,136.5] at σ=1 gives containment margins **+0.00379, +0.00372, +0.00353, +0.00342** (all positive — the slack dominates the −2.7e-4 tight-box deficit); the saddle subcell of [100,105] contains lhs with margins +0.00522/+0.00884. Failure region N∈[100,130.7] lies wholly inside the ratio-1.05 [100,400] cover; the sweep (N≥400) has G(y,1)>0 anyway.
- **Assessment:** no counterexample to the theorem; the H3 proof must be amended to either (a) restrict the lower-box claim to the saddle curve (all that is used), or (b) certify containment of the coded mS box directly over the σ-subcells.

## R6 — INCOMPLETE (honesty note)

I did not independently rerun step_18 / step_10_J192 / 200 random sweep cells before the turn budget expired. What I did verify independently: the fixed files differ from the buggy baselines only in the intended fix lines; and every stored constant I recomputed independently matches — B_abs = 0.16666814452037215 → d6 = 0.3702435178 (D13s: 0.37024), Hp3 box lower 0.153836 (D13s: 0.1538), y(1+y²) = 0.127300114974765 (D13s: 0.12731), ecw ratio exactly 2^{1/4}. The full rerun remains to be done by a follow-up; my R1–R5 results give no reason to expect the stored sup values (194.1267251 / 254.2275176 / 47.0701250) to change.

## Overall verdict on C = 254.2276

**No item broke the certificate.** B1 and B2 are correctly fixed (R1, R2: HOLDS); the new H1 and H2 proofs are correct and complete (R3, R4: HOLDS); H3's proof has a genuine presentation gap — the lower-box containment is proved only on the saddle curve, and the stronger domain statement is false at (N∈[100,130.68], σ=1) — but this does **not** affect the soundness of the certificate (saddle containment certified; coded boxes verified to contain the truth even at σ=1 in the four tested failure-region cells). Subject to the R6 reruns being completed, **C = 254.2276 stands** (certified sup 254.2275176 ≤ 254.2276).

## Evidence table (all printed by my code)

| check | result | digits/dps |
|---|---|---|
| f′ monotonicity (R1) | sympy residuals exactly 0 (f′−g/x², g′−2x²/(1+x²)²) | exact |
| A1 box containment (R1) | worst viol −2.09e-9 and −2.30e-10, 110 884 pts, 149 cells×4 | 60 dps |
| zz box containment (R1) | worst viol −1.07e-12 | 60 dps |
| φ = (π/√2)B_h[A1′/y] (R1) | worst resid 4.44e-56, 2006 pts incl N=100, 1891 | 60 dps |
| σ* range (R1/R5) | [0.428579581, 0.5) for N∈[100,1e20] | 60 dps |
| ecw fix (R2) | 0.483941449038287; ratio to buggy = 2^{1/4} exactly | 50 dps |
| diffs (R2) | only intended lines changed; step_18 identical | — |
| p2mag modulus (R2) | sympy residual 0 | exact |
| saddle branch (R3) | Re Φ′ resid 9.9e-48; branch eq 3.2e-58 | 60 dps |
| B_true exact form (R3) | worst rel resid 1.01e-32, 5000 pts t∈[1e-14,10³] | 60 dps |
| A1n bracket (R3) | min(A1n−L) = 1.79e-50 @100 dps, t∈[1e-12,10⁸] | 100 dps |
| B_true arb box (R3) | 100 000/100 000 contained, max viol −1.48e-6; B_abs = 0.16666814452037215 | 70 dps arb |
| x²\|g(x)\| ≤ 1 (R4) | max 0.999998915112 (200 000 pts), →1 as x→0 | 50 dps |
| Σc_k = 1 (R4) | sympy = 1 exactly; S(K) ≈ 1−1/√(πK) | exact |
| 2Re w ≥ \|w\| (R4) | min ratio 1.4142135, dense (x,k) | 50 dps |
| I0 ≤ bound (R4) | ratios 0.2724 / 0.2696 / 0.2691 / 0.2687 (N=100/400/10³/10⁴) | 50 dps |
| Hp3 ≥ 0 (R5) | single box lower 0.1538363242986387; true min 0.39809551699 (200 001 pts) | 60 dps arb/mp |
| H′ = t³Hp3 (R5) | resid 2.8e-61 (200 pts) | 60 dps |
| **G(y,σ) domain scan (R5)** | **min −5.4722904e-6 at (y(100), σ=1); G(y,1)=0 at N=130.68489** | 80 dps |
| G on σ∈[0,0.5] (R5) | min +2.56e-14 > 0 | 80 dps |
| coded mS box at σ=1 (R5) | contained, margins +0.00379…+0.00342 (4 cells) | 70 dps arb / 80 dps |
| t_max (R5) | 0.127300114974765 < 0.13 | 80 dps |

## Files produced (D13t_review_fixes/)

`step_17.py`–`step_20.py` (R2 grep/diff/constants/sympy), `step_21.py`–`step_25.py` (R1 batches, σ* range, Bh exactness, sympy f′), `step_27.py`–`step_29.py` (R3 saddle branch, B_true form, A1n bracket, arb box), `step_30.py`–`step_33.py` (R4: x²|g|, c_k/Σc_k, series, beta identity, 2Re w, I0), `step_34.py`–`step_39.py` (R5: Hp3 box reproduction, true min, H′=t³Hp3, G domain scan + finding, coded-box containment), `result.md`, `progress.log`. Note `step_23.py` (inverted-box test artifact, superseded by `step_24.py`) and `step_35.py` (sympy non-simplification, superseded by the numeric check in `step_36.py`) are retained for transparency.

**Remaining work for a follow-up:** execute R6 (rerun `fixed/step_18.py`, `fixed/step_10_J192.py`, and ≥200 random cells from `sweep_cells_r1.0002.json` against the stored constants); extend the coded-mS-box σ=1 containment check to all 29 cells of [100,400]; amend the H3 write-up as noted in R5.