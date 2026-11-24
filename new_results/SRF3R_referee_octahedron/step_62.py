report = r"""# SRF3R referee report — SRF3 octahedron end (claims A–D)

Referee: SRF3R (adversarial). Subject: `research/08_theory/SRF3_octahedron_radius_proof/{proof.md, result.md, verify_srf3.py, SRF3_hi_coeffs_K120.json}` and the independent `SRF1_radii_3d/SRF1_n3_K60.json`. All numbers below were printed by SRF3R's own code (`step_21–61.py`, consolidated in `verify_srf3r.py`, 136 s, all checks pass). Notation: z0 = (−log2 + iπ)/2, R = |z0| = 1.6085752558559049, θ = arg z0 = 1.7879528081463, b_k = a_k/k!.

## Verdicts

| claim | verdict | one-line reason |
|---|---|---|
| A (exact representation + 3 forms) | **PROVED** | algebra + Jacobians exact; my direct quadrature (graph & radial, not via T) agrees 50–67 digits at d = ±0.01, ±0.03 |
| B (Borel singularities exactly z0, z̄0, square-root) | **PROVED AFTER REPAIRS** | §3.2 contour claim false for moderate x (V-split repair written); zero/corner/Hankel content independently verified |
| C (radius zero, limsup = 2/√(π²+log²2)) | **PROVED AFTER REPAIRS** | §3.2 GAP repaired by V-split; §3.3 R2-split defective (H0 off by e^w; Ψ-part keeps z0-strength) — replaced by direct non-cancellation; §3.4 correct |
| D (growth law, both series) | **PROVED AFTER REPAIRS** | τ-law and 2nd d-display correct (confirmed 9–10 digits); 1st d-display sign WRONG (planner finding reproduced); fits: β = −1.500 sharp, Cτ 9 digits, Cd 9.7 digits |

Claim E (Borel summability, NUMERICAL): independently confirmed and strengthened — my Padé–Borel [55/55] gives 46.0, 47.3, 46.5, 34.7 digits at d = 0.01, −0.01, 0.03, 0.1 vs my own quadrature (SRF3 reported 37.0, 15.2, 6.3 with [29/30]-type orders).

## Task 2 — sign convention RESOLVED
SRF1's `tau` and `d` are **bit-identical** to SRF3's `eta` and `d` for k ≤ 60 (same variable, same series; step_21). Raw ratios (step_21, verify_srf3r R3b):
- b_k(d) / (k^{−3/2} Im[z0^{1/2−k}]) for k = 52–60: −1.525, −1.599, −1.537, −1.128, −1.547, −1.500, −1.559, −1.525, −1.585 — oscillating about **−√(3π/4) = −1.53499** (the planner's −1.54 reproduced). With Im[z̄0^{1/2−k}] the ratios are +1.52…+1.60.
- b_k(τ) / (k^{−3/2} Re[z0^{1/2−k}]) oscillates about −√(3π/8) = −1.08540 as claimed.
Resolution: SRF1's d IS the same variable; the mixing factor e^{−z̄0} = √2 i is correct; **Claim D's first d-display has a sign error**: the correct statement is
  b_k(d) ~ −√(3π/4) k^{−3/2} Im[z0^{1/2−k}] = +√(3π/4) k^{−3/2} Im[z̄0^{1/2−k}] = +√(3πR/4) R^{−k} k^{−3/2} sin((k−½)θ).
proof.md's **second** d-display (the sin((k−½)θ) form) is correct and confirmed to 9.7 digits (my fit C2 = 1.94682206573 vs theory √(3πR/4) = 1.94682206555; C1 → 0).

## Line-by-line ledger (proof.md)
- §1 Thm A: **CORRECT** (all exponents, the p^{−2} = (1−τ)² prefactor; square/Laplace/logit forms: algebra checked, numerics 24–45+ digits my code; 40 digits theirs).
- §2 Borel contour B_f = (1/2πi)∮f e^{z/τ}dτ/τ: **CORRECT** (residue check).
- §2 Hankel B_f ~ −iE(τ*)²/(2√π) z^{−3/2}e^{z/τ*}: **CORRECT up to the square-root determination** (re-derived exactly via ₁F₁(−1/2;1;w): B ~ −(−τ*)^{1/2}(τ*)^{3/2}/(2√π) z^{−3/2}e^{z/τ*}).
- §2 zeros τ* ~ z0/x, ladder, s=1/2 minimal: **CORRECT** (my expansion: τ* = (z0/x)(1 − log(s(1−s))/(2x) + O(x^{−2})) — §4.1's display has "+" instead of "−" before log(s(1−s))/(2x): **MINOR FIX (sign typo)**; the downstream 1/τ* formula and the [s(1−s)]^{z/(2z0)} factor use the correct sign).
- §3.1 interior: **MINOR FIX** — the integrand is *not* entire in τ even for t_i ≥ δ (Σt_i^{2τ} vanishes, e.g. 2t1^{2τ} = −t3^{2τ} at |τ| = R/log(t1/t3)); but its zeros are uniformly at |τ| ≥ R/log((1−2δ)/δ) > 0 (same fiber equation), so the piece's coefficients grow geometrically and its Borel contribution is entire. Conclusion correct.
- §3.1 corners: fiber equation e^{−2zw}+e^{−2z(1−w)} = −1 correct; min |z| ≥ R with equality only at w→0,1; = 2R at w = ½. **Verified independently** (all-branch search, 13 w-values: min 1.0076R at w = 0.005, 2.0 at w = ½; actual-Q corner scan: min |τ*|U/R ≥ 1 for U ≳ 8). Full uniformity write-out remains **proof-sketch+** (admitted open item; double roots only at w = 1, so roots are simple/continuous on compact w-sets — repair path in FINAL).
- §3.1 kernel |k_j(V)| ≤ (C(V)e^{−V})^j/j!, C(V) = 2e^V|log(1−e^{−V})|: **EXACT** (k_j = L^j/j!); C(V) decreasing 4.17 → 2 verified.
- §3.2 contour "inside all zeros of Q": **GAP (FALSE as written)** — min_{s∈[δ,1−δ]} |τ*|x/R dips to 0.72–0.85 for x ∈ [1.2, 3] (my scans; at s ≈ 0.8–0.95, inside the strip), and to 0.61 at s = 0.015. SRF3's step_06 sampled only s ∈ {0.5, 0.25, 0.1} and at u = 0.1 hit a wrong Newton branch (printed 0.0056), so the failure was invisible to them. **Repair (written out in FINAL)**: split the strip at V0(ε) — for x ≥ V0 all zeros are outside |τ| = (R−ε)/x (verified: min ratio 1.16 at x = 5, 1.036 at x = 20, min at s = ½, asymptotics 1 + log2/x); the bounded-V part contributes an entire function to B_T (uniform zero-free disc by compactness, Q(0) = 3). Same split for the corner pieces with U ≥ U0. Conclusion (limsup ≤ 1/R) intact.
- §3.2 |Q| ≤ 3u^{−2r} chain and the Re τ = +r bound u^{−r}4^{2r}: **MINOR FIX** — the chain needs s(1−s) ≥ u, false in general; correct bounds: |f| ≤ √2[u(1−u)]^r + [(1−u)²s(1−s)]^r ≤ √2+1 on Re τ = −r, and (δ(1−δ))^{−r} replaces 4^r on Re τ = +r; max over the circle at σ = ±r is valid by log-convexity. Constants only; estimate |a_k| ≤ C'k!/(R−ε)^k stands.
- §3.2 "(corner/kernel) same bound with room to spare (2R instead of R)": imprecise — corners attain R/U at edge-ward directions (2R only at the fiber midpoint); harmless for the upper bound. **MINOR**.
- §3.3 H0(w) = e^w√(1+2e^{2w}) as "strict V→∞ limit": **WRONG** — verified numerically Φ(w/V, V) → √(1+2e^{2w}) = e^w√(2+e^{−2w}) (ratio Φ/H0 → e^{−w}: 0.7499 at w = 0.3, 0.5005 at w = 0.7).
- §3.3 "Ψ-part analytic in |z| < 2R": **FALSE as stated** — the Ψ-part keeps z0-strength (πi/8 − 1)·i√2·z0^{1/2} ≠ 0: per-edge strength (πi/4)e^{z0}z0^{1/2} vs Ĥ-strength i√2·z0^{1/2}, ratio exactly πi/8 (verified numerically to 40 digits). Cause: the s-integral's zero shift produces the B(3/2,3/2) = π/8 factor; freezing s-factors at τ = 0 loses it. **Non-cancellation conclusion nonetheless PROVED** by the direct argument (§4.1 computes the total D ≠ 0; bounded-V part entire; higher Hankel terms ∝ (1−z/z0)^{3/2+} cannot cancel the (1−z/z0)^{1/2}; next ladder zeros at |z| = ½√(log²2+9π²) = 4.7253 > 2R). Ledger "K: PROVED (R2 mechanism)" is overclaimed: the |k_j| bound is proved, the corollary is not (as stated); repair = the direct argument. Combinatorial steps (n−j)!/n! ≤ 1/j! (⟺ C(n,j) ≥ 1) and (j+1)^{n−j+1} ≥ 2^{n−j+1} verified.
- §3.4: **CORRECT** — 3-consecutive-k lemma verified to k = 2·10⁵ (weakest window max|cos| = 0.7796 ≥ sin((π−θ)/2) = 0.6263; no irrationality needed, good-gap θ > bad-gap π−θ); |a_k/k!|^{1/k}R ≈ 0.93 at k = 120 with the expected −(3/2)(log k)/k approach. τ↔d biholomorphy: mixing identity verified to 2e−69 (all k ≤ 120), mixing factor e^{−z̄0} = √2 i ≠ 0 preserves the limsup.
- §4.1 singular-part computation: **CORRECT** — E → 2e^{z0}√x (Q'(τ*) = 4x+O(1) since u^{2τ*} → −2; e^{τ*M} → e^{z0}); e^{z/τ*} = e^{zx/z0}[s(1−s)]^{z/(2z0)}(1+o(1)); s-integral → B(1+z/(2z0),·) → B(3/2,3/2) = π/8; x-integral → Γ(−½)(1−z/z0)^{1/2} = −2√π(1−z/z0)^{1/2}; per edge (πi/4)e^{z0}z0^{1/2}(1−z/z0)^{1/2}; ×3 edges; S/S(1) = (2/√3)(1−τ)²T with Borel[τg] = ∫Borel[g] (primitives raise branch order) ⟹ strength (√3πi/2)e^{z0} = −π√(3/8) in the tracked determination; the real form +π√(3/8)[z0^{1/2}(1−z/z0)^{1/2} + c.c.] is fixed by the fit (admitted open phase step) and yields the confirmed law.
- §4.2 constants: **CORRECT** — A = 1/R = 0.6216681478591508, β = −3/2, C = √(3πR/8) = (3π)^{1/2}(π²+log²2)^{1/4}/4 = 1.3766110843131, C_d = √2 C.
- §4.3 invariance: **CORRECT** — T(0) = √3/2, 1/T(0) = 2/√3; my w-series fit: C_w/C_τ = 1.00000000 (8 digits); w₂…w₆ = paper thm:oct values to 2e−50; d²d³d⁴ exact to 1e−71.
- §4.4 table: reproduced independently — β = −1.500 sharp (resid 2.6e−14 at −1.5 vs 1e−8 at ±0.002); Cτ = −1.3766110829 (err 1.4e−9, k0 = 100, 5 corr.); Cd = 1.94682206573 (err 1.8e−10); R, θ profiles sharp at the theory values (resid ×10⁵ at ±1e−5 rel.); vertex-angle amplitude ratio ~5e−13 (their bound <2e−3 confirmed); "φτ = π−θ/2 = −2.2476163" is a sign-convention typo (π−θ/2 = +2.2476162; the fit value −2.2476 is −(π−θ/2), the cos(kθ−φ)-convention); "φd − φτ = 3π/2" is convention-consistent (equals +π/2 in the cos(kθ+φ) convention; both give the confirmed sin((k−½)θ) law).
- §5: NUMERICAL — independently confirmed and strengthened (see above).

## Evidence table (SRF3R's own prints)
| check | result | digits/values | step |
|---|---|---|---|
| S(1+d) graph & radial vs 8(1−τ)²T, d = ±0.01,±0.03 | 50.5–67.5 digits | worst 50.5 (verify_srf3r), 55–67 (dps 55) | 23, 58, R1 |
| T equivalent forms vs simplex | 24–45+ digits | square/Laplace exact at 45 dps; logit 24 | 24 |
| zeros tail min_s |τ*|x/R (s∈[.05,.95]) | ≥ 1 | 1.163 (x=5) … 1.036 (x=20) | 27, R2a |
| zeros moderate-x min |τ*|x/R | < 1 (§3.2 GAP) | 0.850 (x=2), 0.72 (x≈1.2–1.7) | 26, R2b |
| corner fiber min |z|/R | ≥ 1; 2 at w=½ | 1.0076 (w=.005)…2.0 (w=½) | 29, R2c |
| actual-Q corner min |τ*|U/R | → ≥ 1 | 0.81 (U=3), 1.006 (U=8), 1.044 (U=12) | 31 |
| d-sign ratio k=52–60 | −√(3π/4) | −1.525…−1.585 vs −1.53499 | 21, R3b |
| τ C1 fit (k0=100, 5 corr) | −1.3766110829 vs −1.3766110843131 | err 1.4e−9 (9 digits) | 39 |
| d C2 fit | 1.94682206573 vs 1.94682206555 | err 1.8e−10 (9.7 digits) | 39 |
| β profile | sharp min at −1.5 | resid 2.6e−14 vs 1e−8 at ±0.002 | 40 |
| R, θ profiles | sharp at theory | resid ×10⁵ at ±1e−5 rel | 38, 40, R4b |
| mixing identity | worst rel err | 2e−69 | 48, R5a |
| d²,d³,d⁴ vs paper | abs diff | ≤ 1e−71 | 42, R5b |
| w₂..w₆ vs paper; C_w/C_τ | abs diff; ratio | ≤ 2e−50; 0.9999999978 | 42 |
| vertex amplitude ratio | < 5e−13 | (claim <2e−3) | 41 |
| 3-consecutive-k lemma | 0 violations | k ≤ 2·10⁵, weakest 0.7796 ≥ 0.6263 | 44, R5c |
| Φ(w/V,V) limit | √(1+2e^{2w}), not e^w√(1+2e^{2w}) | ratios e^{−w} to SRF3's H0 | 33 |
| strength ratio edge/Ĥ | exactly πi/8 | 0.39269908 i to 40 digits | 33 |
| Padé–Borel [55/55] vs own quadrature | 34.7–47.3 digits | d = 0.01, −0.01, 0.03, 0.1 | 43, R6 |

## Remaining open items (both admitted by SRF3; status unchanged: proof-sketch+)
1. §3.1 corner uniformity (smooth partition of unity gluing edge strips to corners; Rouché estimates on compact w-sets; double roots only at w = 1 so the roots are simple — repair sketched in FINAL).
2. Hankel transfer phases (uniform error of the per-(s,u) Hankel asymptotic in the tail; analytic tracking of the square-root determination — currently the sign of D is fixed by the 1e−8 fit).
"""
open("result.md","w").write(report)

prog = open("progress.log","a")
prog.write("[SRF3R] zeros scan: min_s |tau*|x/R = 0.61..0.85 for x in [1.2,3] INSIDE the strip => S3.2 contour claim FALSE as written; tail (x>=5) >= 1.03 at s=1/2. V-split repair works (bounded-V part entire).\n")
prog.write("[SRF3R] corner fiber eq: min|z| >= R all w (1.0076R at w=0.005; 2R at w=1/2), all-branch search; actual-Q corner: >=1 for U>=8. C(V) decreasing 4.17->2, |k_j| bound exact.\n")
prog.write("[SRF3R] S3.3 defective: Phi(w/V,V) -> sqrt(1+2e^{2w}) not e^w sqrt(1+2e^{2w}) (ratios 0.7499/0.5005 = e^{-w}); Psi-part keeps z0-strength (ratio edge/H = pi i/8 exactly). Non-cancellation survives via direct D != 0.\n")
prog.write("[SRF3R] fits (own method): beta=-1.500 sharp; C_tau err 1.4e-9; C_d err 1.8e-10; R,theta sharp; vertex amplitude 5e-13; mixing 2e-69; w2..w6,d2..d4 exact; 3-consec-k lemma to 2e5 OK.\n")
prog.write("[SRF3R] Theorem A re-verified: graph+radial direct quadrature vs 8(1-tau)^2 T: 50.5-67.5 digits at d=+-0.01,+-0.03 (exact decimal p; float-1.01 pitfall identified, dS/dp~6.9 explains 6.1e-17).\n")
prog.write("[SRF3R] Pade-Borel [55/55] vs own quadrature: 46.0/47.3/46.5/34.7 digits at d=0.01/-0.01/0.03/0.1 (stronger than SRF3's table).\n")
prog.write("[SRF3R] verify_srf3r.py written, all PASS in 136s. VERDICTS: A PROVED; B,C,D PROVED AFTER REPAIRS (repairs written in result.md/FINAL); claim D first d-display sign error CONFIRMED (planner's -1.54 reproduced).\n")
prog.close()
print("done")
