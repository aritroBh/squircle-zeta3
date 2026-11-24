
rmd = """# DEC5D — the regular decagon's P⁽⁵⁾(1): depth two needed? — status report

## Bottom line (preliminary)

Value computed to ≥113 digits (two quadrature levels); exact edge integrand derived and
assembled (matches to 101 digits); letters {0, 1, −φ, φ²}; depth test in progress.

## Task 1 — NUMERICAL: P⁽⁵⁾(1) and sanity checks (done)

Engine `dec_engine.py` (Route A d-series of polygon_engine.py, one wedge (−π/10,π/10), ×10 symmetry):
- two independent runs dps 120 / dps 170 agree: |ΔP⁽⁵⁾| = 1.6e-114 (P⁽⁰⁾..P⁽⁴⁾ agree to 116–126 digits);
  final tanh-sinh level diffs 1.08e-117 (dps 120), 1.85e-167 (dps 170).
- P(1) = 2.008114158862272798697976726327421423305944607218666812… = 10L, L = √((25−11√5)/10): residual 2.9e-177.
- P‴(1) = −4.776200180882238416632022453444885836970773916879016004… = 60·L·S with the paper/DEC1
  15-term ℚ(√5) form S: residual 6.6e-170.
- **P⁽⁵⁾(1) = −10.1644521920838037189075927518683458765038789365773310500737349139543150766489835168246984814886990478204766055284955…**
- P⁽⁵⁾(1)/L = −50.616904159684958298011819837155762859466441700415750720050351864156418026604209915752195581787514539343384389364197798…

## Task 2 — exact structure (done)

Edge (vertical, x = −h): V0 = (−h, R/2φ), V1 = (−h, −R/2φ), h = 1/(2φ), R = 1/(2φ² sin(π/5)).
`symeng.edge_params` (numeric) gives j0 = 2, j1 = 3, nonloc = {0, 1, 4}, α = β = 1/(2φ²), L, s0 = −L/2.
Exact ℚ(√5) data (derived, matches all 50-digit values):
- local: a₂ = t/(2φ²), a₃ = (1−t)/(2φ²);
- non-local: a₁ = (1+t/φ)/(2φ²) (κ₁ = 1/φ), a₄ = (1−t/φ²)/(2φ) (κ₄ = −1/φ²), a₀ = 1/(2φ) CONSTANT;
- s0/L = −1/2, h²/L² = (5+2√5)/4; E5 = L·(polynomial with ℚ(√5) coefficients).
Exact substituted E5/L vs `symeng.E5_sym` (conv='B'): max diff 2.3e-50 at 5 t-values (dps 50).
Normalization P⁽⁵⁾(1) = 5!·10·∫₀¹E5 = 1200∫E5: residual 2.7e-41 (dps-40 quad).
Expansion: 461 monomials (a,b,c,d,e,f) in (log t, log(1−t), w₁=log(1+t/φ), w₄=log(1−t/φ²), log2, logφ),
degree counts {1:6, 2:21, 3:56, 4:126, 5:252} — the octagon's exact structure; all coefficients L·(p+q√5)·tⁿ, n ≤ 5.
Assembly from 756 numerically integrated monomials (dps 110): matches engine P⁽⁵⁾/L to **1.6e-101**.
Weight-5 part W5/L = −450.306606005917383005102743374749276889885908477643451891… blocks:
local −495.637734500553…, single-1 −93.285053269939…, single-4 +109.034320420277…, mixed +29.581861344298… .

## Task 3/4 — letters, V, depth test: IN PROGRESS
"""
open("result.md","w").write(rmd)
print("result.md updated")
