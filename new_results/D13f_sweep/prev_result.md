# D13e — Certified decay theorem: W2/W3 fixed and certified at N = 100, 400; all-N sweep architecture built, one identified gap remaining

## Bottom line

**(1) The W2/W3 `nan` bug is FIXED and both N=100 and N=400 now carry fully arb-certified assemblies** (all 9 terms + sum C, printed as intervals). The certified constants are **C(100) = 42.425, C(400) = 46.545**, better than D13d's estimate of 60–75. The certified bounds dominate the mpmath "truths" for every piece.

**(2) The all-N≥N₀ sweep is architected but NOT finished** (one interval-arithmetic gap precisely identified, see §3). So the statement "for all N ≥ N₀" is **not yet a theorem**; the single-N statements at N=100 and N=400 **are** certified (each inequality below is an arb enclosure).

## 1. The fix (task item 1) — PROVED for the two evaluated N

`step_34.py` (new, this session). Three repairs over step_33:

- **γ via lgamma difference** (`logg`): γ(x)=log g(x)−log g(x₀) with `log g(x) = log(√π/4) + lgamma(s/2) − lgamma(s/2+3/2) − log x`, s = e^{iπ/4}x. Since Re(s/2) = x/(2√2) > 0, this is the analytic branch — no cut crossing. Verified against 60-digit mpmath: contains truth.
- **Midpoint + certified reach** on each of J=24 subintervals (kills interval dependency): |R(u)| ≤ |R(uₘ)| + M₃c·u₁²h/2, |R₄| ≤ |R₄(uₘ)| + M₄c·u₁³h/6, |γ| ≤ |γ(uₘ)| + G₁c·h, |γ₂| ≤ |γ₂(uₘ)| + G₂c·u₁h, with M₃c,M₄c from the regular explicit Φ‴/Φ⁗ bounds and G₁c,G₂c from the elementary bounds ψ′(t) ≤ 1/t²+1/t, |ψ″(t)| ≤ 2/t³+1/t² (integral test).
- **Curvature minima by cell coverings** (rigor fix: step_33 used point grids, which do not bound the min over an interval).
- **W3 q-factor**: |e^z−1−z| ≤ |z|²·max(1,e^{Re z})/2 with e^{Re(R+γ)} ≤ e^{Re R}|g(x₀+u)/g₀| bounded per cell via |Γ(x+iy)| ≤ Γ(x) (elementary from the integral definition).

Certified output (arb intervals, ctx.dps=70; `/sc`, sc = N^{−3/2}e^{−√(πN)}):

| N | W1 | W2 | W3 | TM | TR | TLm | TLf | TL0 | D | **C** |
|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 10.7454 | 1.8813 | 11.0078 | 11.2800 | 7.9743 | 17.8193 | 5.3823 | 1.5e−18 | 0.3503 | **42.4247** |
| 400 | 11.9079 | 4.4679 | 28.0595 | 3.9607 | 8.9461 | 14.8142 | 0.4559 | 1.5e−90 | 0.3185 | **46.5448** |

(Full 60-digit intervals in `cert_single_N.json`.) Each number is an arb enclosure; the C printed is `((2/π)(W1+…+TL0)+D)/sc`, certifying
|M_k^main + (√π/N)e^{−√(πN)}cos√(πN)| ≤ C·N^{−3/2}e^{−√(πN)} at N=100 (C=42.43) and N=400 (C=46.55).

**Cross-check vs truth** (`step_35.py`, mpmath 60dps quadrature of the exact pieces): W1 true 2.03/3.39, W2 true 0.56/0.69, W3 true 1.49/1.13 (N=100/400) — every certified bound dominates its true piece.

## 2. All-N sweep (task item 2) — architecture built, one gap left

`step_36.py` implements the box assembly with the machinery needed for wide N-boxes, all derived and numerically verified this session:

- **σ-parametrized saddle enclosure** ρ := y₀/y* = 1+y*²σ, σ∈[0,1] — this *is* D13d's proved enclosure y₀∈[y*, y*(1+y*²)]; I re-proved the two polynomial inequalities self-containedly (F(t)=t²(1+√2t)/(1+√2t+t²) increasing since 2+4√2t+4t²+√2t³>0; F(y(1+y²))≥y² ⟺ 1+√2t(2+y²) ≥ y²+y⁴, true for y≤0.78 i.e. N≥2.6) — the correlation ρ=1+y²σ is restored exactly where interval arithmetic loses it (ω's numerator, the D-phase B̂ = B/y).
- **Cancellation-free regular forms**, all derived and numerically verified: curvature |ReΦ″|/2 = (πρ²y/4)·(√2η³ρ+4η²ρ²y+2√2ηρ³y²+ρ⁴y³)/(η²(η²+√2ηρy+ρ²y²)²) at η=x/x₀ (identity A₂=yB₂ re-verified); p₂, 1/p₂, |p₂| regular; e^{ReΦ(x₀)+θ} = e^{−(N/2)S(y₀)} with −S(t)/t³∈(√2/3)[1/(1+√2T+T²),1]; TLf/TL0 in e^{−θ}-canceled regular form; D-phase φ = ImΦ₀−θ = (π/√2)B̂(y,σ) regular.
- Point-N boxes reproduce C ≈ 43 (N=100), 46.6 (N=400), consistent with step_34.

**The one remaining gap (precisely identified):** over the widened (y,σ)-boxes, the *complex* digamma/lgamma evaluations (`gapF`, `logg`) blow up by cancellation (|γ′(x₀)|-box came out ±41 instead of ≈0.13, inflating W1), and `λR = −ReΦ′(xr)` loses its sign over wide boxes. Fix (same as step_30's and as used for G₁c/G₂c): for box evaluations replace all γ-point-evaluations by the **elementary real bounds** |γ′(x)| ≤ (3/4)ψ′(x/(2√2))+1/x ≤ (3/4)(1/t²+1/t)+1/x etc., and compute λR via the regular G-form λR = a−(πρ/2)Ĝ(η) (derived this session). Then re-run the sweep: cells [100,400] (ratio ~1.5), dyadic to N_J≈10¹⁰, final cell y*∈[0,y_J] with the regular final-cell assembly (h₁ two-sided via the Weierstrass product + Euler–Maclaurin bounds I derived and numerically validated: h₁(y₀)→1 with h₁−1≈−0.53y₀).

## 3. Status vs the task list

| item | status |
|---|---|
| (1) Fix W2/W3; N=100,400 with all 9 terms + C as arb intervals | **DONE, certified** |
| (2) Dyadic N-boxes + analytic tail | architecture + regular forms done; sweep blocked on the digamma-over-boxes gap (fix identified); **not certified** |
| (3) Table check vs mk_deg9.tsv | not run (bound ≈ 42–47·sc vs truth ≈ 0.6–1.6·sc from D13d, consistent) |
| (4) Rest term | not started |

## 4. Evidence table

| check | result | residual | digits |
|---|---|---|---|
| lgamma-difference γ vs mpmath (N=100 cell) | contains truth | inside [±0.62] box pre-subdivision | 60 |
| step_34 C(100) | [42.42472910… ± 8.1e−64] | arb width 2e−62 rel | 62 |
| step_34 C(400) | [46.54477870… ± 5.6e−63] | arb width 1e−62 rel | 62 |
| W1,W2,W3 truths (quad) vs certified, N=100 | 2.03/0.56/1.49 ≤ 10.75/1.88/11.01 | domination | 60 |
| W1,W2,W3 truths vs certified, N=400 | 3.39/0.69/1.13 ≤ 11.91/4.47/28.06 | domination | 60 |
| Binet ψ-bound, Γ-ratio bound, h3 enclosure | margins ≥ 0 on 200-pt grid | <1e−37 | 50 |
| h₁(y₀)→1, product formula, diagonal Stirling rate | match | <1e−60 | 60 |
| curv2 η-form vs direct ReP2 | match (0.1375 vs 0.1384) | interval overlap | 70 |

## 5. Files produced this session

- `step_34.py` — **the certified single-N assembly (main result)**; output in `cert_single_N.json`.
- `step_35.py` — mpmath cross-check: certified W1/W2/W3 dominate exact pieces.
- `step_36.py` — y-form box assembly (σ-enclosure, regular curvatures/phase; has the identified digamma-over-boxes gap for wide boxes; correct on narrow boxes); `cert_boxes.json` (partial).
- `dbg1.py` — diagnostic isolating the gapF box blow-up.

**Estimate for the record (unchanged):** with the box sweep completed as designed, C(N) ∈ [42, 55] for all N ≥ 100, against a true constant ≈ 2. What is PROVED today: the inequality at N=100 with C=42.43 and at N=400 with C=46.55, every step an arb enclosure.