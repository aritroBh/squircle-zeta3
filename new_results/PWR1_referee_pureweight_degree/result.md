PWR1 referee report — Theorem thm:pw (pure weight) and Proposition prop:degree (degree of the Fourier coefficients)

## Scope and method

I re-derived every identity in both written proofs independently, rebuilt all exact data with my own code (formal-zeta engines `pwr1_engine.py`, `pwr1_rhombus.py`, hand-rolled rational series for the Gudermannian table), and ran fresh quadratures at 50–105 digits. Full report: `result.md`; verifier `verify_pwr1.py` (106 s < 5 min, **21/21 checks PASS**, output `verify_pwr1.out`).

## Verdicts

- **thm:pw: PROVED.** Every step of D16's mechanism (M1–M3, single-appearance lemma, Rat_n, Λ_n, the n-dimensional part via H1's T6, the rhombus case) is correct; my independent engine reproduces D16's, H1's, and D11's exact data with difference exactly 0; γ cancels exactly at every order; every coefficient is weight-homogeneous; all interchanges are justified at the level of formal Taylor coefficients.
- **prop:degree: PROVED.** K3's master identity, layer identity, Watson/D̃ formula, degree bound, and both top-degree closed forms are correct; the degree is attained (tops manifestly nonzero); the leading-coefficient formula matches the paper; the "uniformity in k" concern resolves cleanly because the m-sum at each fixed order τᴺ is finite.

## Repairs (documentation level only; no mathematical error found)

- **R1 (MINOR FIX, K3 result.md):** the displayed formula I_{m,j} = j!m!Σ_{l≥j}(−1)^{l−j}s(l,j)/(l!·l^{m+1}) is correct only with *unsigned* Stirling numbers of the first kind; with signed s(l,j) the factor (−1)^{l−j} must be dropped. (Their code used quadrature, so no numerical consequence.)
- **R2 (MINOR FIX, D16 result.md):** "4u_k(2) = κ_k" is stated/checked only for odd k = 3..11; for even k the η-vs-τ sign convention flips the sign (κ₂ = −1 vs 4u₂(2) = +1; κ₄ = −1/8 vs 4u₄(2) = +1/8). State the convention.

## thm:pw — step table

| step | verdict | key printed evidence |
|---|---|---|
| H1 T1 representation ≡ D16's B_n form (Gamma(1+a), a=−(n−1)τ/n) | CORRECT | re-derived (Dirichlet/Lukacs + size-bias shift, 1+a−τ/n = ε); quadrature n=2,3 matches series |
| M1 moment formula + pattern counting | CORRECT | my engine implements only this formula; reproduces D16 `tau_series_2d.json` through τ¹⁰, H1 routeA (n=2,3) and logB through τ⁸, D11's f_n through η⁸ — all differences exactly 0 |
| M2 (ζ(k) born paired with τᵏ) | CORRECT | weight-homogeneity verified exactly: True at every order k≤10 (n=2), k≤8 (n=3..6) |
| M3 γ-cancellation (moments: Σν_v = s; prefactor: −nγa = γ(n−1)τ vs −γ(n−1)τ) | CORRECT | formal GAM tracked: coefficient exactly 0 at every order, n=2..6 |
| Rat_n; Rat₂ = 1; Rat₃ = (1−τ)/(1−2τ) = 1/(2−p) | CORRECT | re-derived both forms; engine: logRat(2) = 0, logRat(3) = log(1−τ)−log(1−2τ) |
| Λ_n = −Σζ(k)(n−1)ᵏ(1−n¹⁻ᵏ)τᵏ/k | CORRECT | re-derived, matches engine |
| w₁ = 0; B_n\|_{Z=0} = 1 | CORRECT | Σ_s binom(j,s)(−1)^{j−s}n^{j−s}n^s = n^j(1−1)^j = 0 re-derived |
| single-appearance lemma; u_k(2) | CORRECT | u_k(2) = 1/4, 1/4, 1/32, −3/16, 91/192, 129/64, −5039/1024, −22805/768, 561331/5120 (k=2..10), matches D16 exactly |
| n-dimensional part (H1 T6: A₂, A₃, W₄) | CORRECT | e.g. w₂(4) = −21Z₂/8 = A₂(4)ζ(2), w₃(4) = −125Z₃/16 = A₃(4)ζ(3) |
| rhombus (D11c representation + ℚ[q²] purity) | CORRECT | representation verified at 50 digits (4 parameter sets); my q-symbolic engine: pure weight + ℚ[q²] exact through τ⁶; w₂ = (1−q²)ζ(2)/4, w₃ = (1−q²)ζ(3)/4; quadrature fits: w₂ 11 digits, w₃ 7.5–8 digits, w₄ 3.8–5.4, w₅ 2–3.3 |
| paper's displayed series through τ⁵ | CORRECT | matches engine exactly (through τ¹⁰) |
| interchanges (limits/sums/integrals) | CORRECT | everything is finite at each fixed order τᴺ; log-moments of Gamma exist to all orders |

**Numerics (printed):** 100-dps quadrature at 9 τ-values: remainder ratio R/(w₉τ⁹) = 0.99704 → 1 and (R−w₉τ⁹)/(w₁₀τ¹⁰) = 1.00277 → 1 as τ→0 (validates w₂..w₁₀ to ~30 digits collectively); LSQ extraction of w₂..w₁₂: 20.5, 16.8, 13.9, 11.2, 9.0, 6.9, 5.2, 3.7, 2.4 digits (asymptotics-limited); purity residuals w₂−ζ(2)/4 = −1.17e-21, w₄+19π⁴/2880 = −8.2e-15, w₅+(ζ2ζ3/2+3ζ5/16) = +7.3e-12 (γ/rational/lower-weight pollution excluded at these levels); n=3 simplex quadrature remainder ratio 1.0016 → 1.

## prop:degree — step table

| step | verdict | key printed evidence |
|---|---|---|
| a_k − A_k = (−1)ᵏB(τ,2k), B even entire in κ | CORRECT | re-derived; A_k = −(−1)ᵏ4√2/(π(16k²−1)) verified to 7.4e-13 (k=0..5) |
| [τᴺ]B(τ,2k) finite m-sum (no uniformity-in-k issue) | CORRECT | D(τ;m) = O(τ^{2m+2}), so m ≤ N/2−1 only |
| layer identity h_p−h_1 = (1+e^{−2y})^{−1/2}(e^{τL(y/τ)}−1) | CORRECT | re-derived exactly |
| Watson expansion; D̃ formula; e_{r,m} ∈ ℚ | CORRECT | my exact rebuild matches all 50 entries of K3's table; e_{r,2r} = 1, e_{r,2r+1} = 1/2 (proved), e_{r,2r+2} = −5/8−r/3 (r≤5) |
| I_{m,j} alternating MZV of weight m+j+1 | CORRECT | re-derived; quadrature = Stirling series; I_{m,1} = m!η(m+2) and I_{0,2} = ζ(3)/4 exact |
| degree ≤ ⌊n/2⌋−1 | CORRECT | j ≥ 1 ⟹ D_{n,m} = 0 for n < 2m+2 |
| even tops D̃_{2n,n−1} = 2(2n−2)!η(2n) > 0 | CORRECT (PROVED) | diffs 0, 0, 2.2e-103, 1.4e-101, 9.1e-100 (n=2..6) |
| odd tops D̃_{2n+1,n−1} = (2n−2)![ζ(2n+1)+2Ση(2j)ζ(2n+1−2j)] > 0 | CORRECT (PROVED) | classical A_{2m} evaluation verified to ≤1.1e-106 (m≤4); closed form matches D̃ to ≤8.9e-103 (n=1..5); positivity manifest |
| leading coefficient (−1)^{m−1}(2^{2m−1}−1)ζ(2m) | CORRECT | identity 2(−4)^{m−1}(1−2^{1−2m}) = (−1)^{m−1}(2^{2m−1}−1) exact (m=1..8); numeric match ≤2.2e-103 (m=1..6) |
| Θ_n(k²), n ≤ 8 | CORRECT | my D̃-build vs K3's Qpolys: max diffs 1.6e-51, 1.1e-51, 4.5e-49, 3.9e-49, 9.2e-48, 4.0e-47, 2.4e-45 (n=2..8, u ∈ {0,1,2,5,7}); paper's Θ₅ matches |
| a_k quadrature extraction (independent) | CORRECT | extracted Θ_n(k²), k=0..6: rel.err n=2: 6.5e-12, n=3: 3.1e-8, n=4: 1.7e-7, n=5: 6.5e-5, n=6: 1.8e-4 (fresh 5-pt run: 8.9e-10/2.9e-6/2.4e-5 for n=2/3/4); degree-overfit coefficient negligible at every n ≤ 6 |

## Files produced (work dir)

- `verify_pwr1.py` + `verify_pwr1.out` — 21/21 PASS, 106 s.
- Referee engines: `pwr1_engine.py` (formal-Z exact S_n engine with γ/weight tracking), `pwr1_rhombus.py` (q-symbolic rhombus engine).
- `result.md` (full report with re-derivations), `progress.log`.
- Data: `w_exact_n2_tau10.json`, `w_exact_n3_tau10.json`, `quad_logP_n2.json`, `quad_logP_n2_hi.json`, `quad_logS_n3.json`, `quad_rhombus.json`, `quad_ak.json`, `e_rm_mine.pkl`, `I_mj_quad.json`, `fit_w2_w12.json`, `theta_extracted.json`, `tops_check.json`; steps `step_01.py`–`step_75.py` (auto-saved) contain every run that produced the numbers above.

**Bottom line: thm:pw — PROVED; prop:degree — PROVED. Both papers' "not yet independently checked" flags can be lifted after applying the two documentation-level repairs R1, R2.**