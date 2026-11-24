# D15 — ζ(3) at every corner? The diamond-end expansion for p-smoothed 2m-gons

**Working directory:** `new_results/D15_polygon_p1/` (scripts `step_12.py`–`step_33.py`, auto-saved by the harness; key ones listed at the end). Method followed R2 (`R2_rhombus_ab/result.md`) and D11 as instructed; LEDGER "Corner angles" section read.

**Setup.** n_i = (cos(iπ/m), sin(iπ/m)), i = 1..m; K_p = {Σ_i |⟨n_i,x⟩|^p ≤ 1}. m = 2 is the squircle (rhombus at p = 1); K_1 is a regular 2m-gon with circumradius tan(π/2m), corners on the lines ⟨n_j,x⟩ = 0, interior angle α = π(m−1)/m, edge midpoints on the bisector rays.

---

## 1. Exact representation (task 1) — PROVED

**Theorem (representation).** With A = π/m, in polar coordinates measured from a corner direction, the boundary of K_p is ρ(θ,p) = S_p(θ)^{−1/p}, S_p(θ) = Σ_{i=0}^{m−1} |sin(θ + iA)|^p (this is just Σ_i|⟨n_i,x⟩|^p = 1 with ⟨n_i,e(θ)⟩ = sin(θ+iA)). Dihedral symmetry of the configuration gives

$$P(p) = 4m\int_0^{A/2}\sqrt{\rho^2+\rho_\theta^2}\,d\theta,\qquad A/2=\frac{\pi}{2m},$$

one corner-arc (corner → edge midpoint) per fundamental domain. On (0, A/2] every sine is positive, so writing p = 1+d,

$$\log\rho = -e\log S,\quad e=\tfrac1p=1-d+d^2-d^3+\cdots,\qquad S=\sum_i s_i\,e^{d\log s_i},\quad s_i=\sin(\theta+iA),$$

gives P(p) = Σ P_k d^k with **explicit 1D coefficient integrals**: P^{(k)}(1) = 4m·k!∫_0^{A/2} D_k dθ, where D_k is a trigonometric expression in s_i, cos(θ+iA), S₀ = Σs_i = sin(θ+(m−1)A/2)/sin(A/2) times products of ≤ k of the logs log s_i, log S₀ (only integrable log singularities at θ = 0; expansion validity by R2's endpoint-power×log domination argument). The full D_0..D_3 integrands are implemented in `step_23/28/29.py` (T_k = Σ s_i log^k s_i, U₁ = T₁/S₀, U₂ = T₂/(2S₀) − T₁²/(2S₀²), U₃ = T₃/(6S₀) − T₁T₂/(2S₀²) + T₁³/(3S₀³), R₁ = logS₀−U₁, R₂ = −logS₀+U₁−U₂, R₃ = logS₀−U₁+U₂−U₃, Q_k = ∂_θR_k with Q₀ = −S₀′/S₀, D = ρ√(1+Q²) expanded to order 3).

**Verification (NUMERICAL, 50–60 digits).** Corner-arc formula vs direct full polar perimeter (per-segment tanh-sinh over [kA,(k+1)A]): residual **0** at 50 dps (m=3: p = 1.1, 1.3; m=3 p=3: 1.1e-50; m=4 p=2: residual 0); low-precision polygonal chord check agrees to 7 digits. Values (corner-arc, 50 dps):

| m | p | P(p) |
|---|---|---|
| 3 | 0.9 | 3.21978303861926628369269120867573794485687 |
| 3 | 1.1 | 3.69840169493827649453206212737602498657709 |
| 3 | 1.3 | 4.12241778509667439523643589746867721410134 |
| 3 | 2 | 5.13019932064745638217614543868402564586904 |
| 3 | 3 | 5.79943965614898293779820587189150802572358 |
| 4 | 0.9 | 2.25470241726970289486780613861477612765760 |
| 4 | 1.1 | 2.80331571619489499530549088257563331573595 |
| 4 | 1.3 | 3.28376785642123185069075139244818329798032 |
| 4 | 2 | 4.44288293815836624701588099006069369861462 |
| 4 | 3 | 5.26727484556007143279082977683159219086809 |

At p = 1 the arc integral evaluates in closed form (S₀ is a single sinusoid, the polar equation of the straight edge): **P(1) = 4m sin²(π/2m)/cos(π/2m) = 4m(sec(π/2m) − cos(π/2m))** — matches the 2m-gon perimeter to 60 digits (m = 2..6; residuals ≤ 6.2e-61). In particular P(1) = 4√2 (m=2), 2√3 (m=3), 16(sec(π/8)−cos(π/8)) (m=4).

## 2. P′(1): a universal first law — PROVED

**Theorem.** Differentiate the representation at d = 0: ∂_d(ρΦ)|₀ = ρ₀(R₁Φ₀ + Q₀Q₁/Φ₀), Φ = √(1+Q²), Q₁ = ∂_θR₁. Integrating the Q₁-term by parts gives
P′(1) = 4m[ρ₀Q₀R₁/Φ₀]₀^{A/2} + 4m∫₀^{A/2} R₁·[ρ₀Φ₀ − ∂_θ(ρ₀Q₀/Φ₀)] dθ.
The bracket simplifies (two lines of algebra) to **S₀(S₀+S₀″)/(S₀²+S₀′²)^{3/2} ≡ 0**, because S₀″ = −S₀ — the polar differential equation of a straight line. Only boundary terms survive; Q₀(A/2) = 0 (bisector symmetry) and ρ₀Q₀/Φ₀|₀ = −sin²(A/2)/cos(A/2), so with V = S₀(0) = Σ_{i=1}^{m−1} sin(iA) = cot(A/2):

$$\boxed{\;P'(1)=P(1)\cdot R_1(0),\qquad R_1(0)=\log V-\frac1V\sum_{i=1}^{m-1}\sin(iA)\log\sin(iA)=\frac1V\sum_i\sin(iA)\log\frac{V}{\sin(iA)}\;}$$

(the by-parts is legitimate: R₁ is continuous at 0, Q₁ ~ −log θ integrable). Since V > 1 ≥ sin(iA) with one strict inequality for m ≥ 3: **P′(1) = 0 ⟺ m = 2** (then V = sin(π/2) = 1), explaining R2's proved P′(1) = 0 and the literature result as the degenerate case. Corollaries: **m=3: R₁(0) = log 2**, so **P′(1) = 2√3·log 2**; **m=4: R₁(0) = log(1+√2) + (1−√2/2)log 2**.

**Verification:** formula vs coefficient integral: residuals 4.2e-81 (m=3, m=4; 80 dps); PSLQ on P′(m=3): relation [1,−2] vs 2√3 log 2, residual 2.5e-90; PSLQ on P′(m=4) over {1,√2,cos(π/8),√2cos(π/8)}×{1,log2,log(1+√2)}: residual 1.6e-90, maxcoeff 80. R₁(0) values: m=2: 0; m=3: 0.69314718055994530941…; m=4: 1.08439169584521473879…; m=5: 1.35816556700394887954…; m=6: 1.56944298142728768188…; m=8: 1.88771788137165964708….

## 3. P″(1), P‴(1) for m = 3 — NUMERICAL (strong)

Coefficient integrals computed at 100 and 170 dps (two independent runs; a typo in the first 170-dps run was caught by the finite-difference cross-check and fixed — the printed values below are the FD-confirmed ones). Cross-checks: m=2 reproduces R2's PROVED P″ = π²√2/3, P‴ = 6√2(ζ(3)−π²/3) to 59/57 digits; m=3 derivatives vs Richardson-extrapolated central differences of the direct perimeter: P′ agrees to 19 digits, P‴ to 11 digits.

**m = 3** (PSLQ at 170 dps, small bases over Q(√3); rejected the trivial Li₃(1/2)-identity red herring and all large-coefficient relations):

$$P''(1)=\frac{\sqrt3}{8}\,\pi^2+2\sqrt3\log^2 2-4\sqrt3\log 2 = -1.0010944397384099296213510588753545017133849299889\ldots$$
PSLQ: basis {1,π,log2,log3,π²,πlog2,πlog3,log²2,log²3,log2log3,Li₂(1/3),Cl₂(π/3)}×{1,√3}, relation (−8, …, +1·√3π², +16·√3log²2, −32·√3log2), **residual 9.9e-90 (100-dps run), maxcoeff 32**.

$$P'''(1)=\frac{\sqrt3}{8}\Big(16\log^3 2+3\pi^2\log 2-96\log^2 2+96\log 2-9\zeta(3)-6\pi^2\Big) = -5.1454550171356158371716811311239750711155535664197\ldots$$
PSLQ: 32-element basis {1,log2,log3,π,π²,log²,loglog,Cl₂(π/3),Li₂(1/3),πlog,ζ(3),π³,π²log,πlog²,log³,Cl·log,πCl,Li₃(1/3),Li₂(1/3)·log,πLi₂(1/3)} (redundant Li₃(1/2) removed) — **residual 1.9e-139 at 170 dps, maxcoeff 96**, relation 8·(P‴/√3) = 96log2 − 6π² − 96log²2 − 9ζ(3) + 3π²log2 + 16log³2.

**Answers to the main questions (m=3 evidence):**
- **P″(1) does *not* stay a pure multiple of π²**: the π²-coefficient is algebraic (√3/8) but log²2 and log2 terms appear — they are forced once P′(1) ≠ 0 (R₁(0) ≠ 0); m=2 was degenerate.
- **P‴(1) DOES contain ζ(3)** with an explicit angle-dependent coefficient: **−9√3/8** (m=3) vs **+6√2** (m=2, PROVED via R2). Per corner: −3√3/16 vs 3√2/2 — the coefficient *changes sign* with the angle. Everything else entering is elementary: π², π²log2, log^k2 (k = 1,2,3). **No Clausen Cl₂(π/3), no Li₂(1/3), no Li₃(1/3), no π³, and log 3 cancels entirely** — the log-sine/Clausen terms drop out exactly as at m=2 (where only ζ-values survive), despite the basis offering them.

**m = 4** (170-digit values computed, identification **open**):
P″(1) = −1.4562051103812046446537306966795936456777152413012283082812637083293482631402138804126453990677125615381393921189383069887917998883431152326…;
P‴(1) = −4.3944963052760521599791101138060999138223774460125653387613594510957926611720826055722134086632599709180883019106072943224285166954307901073…
(saved in `vals_m34_170.json`). PSLQ over {1,log2,L= log(1+√2),π²,log²2,L²,log2·L}×field₄ returned only a 29-term relation with maxcoeff 19354 and residual 2.4e-109 « 1e-170 — **rejected as an artifact**; scalings by c8, √2c8, s8, √2s8, P(1) gave nothing. The m=4 answer needs a wider basis (Cl₂(π/4), Catalan G, Li₂(√2−1), …) and/or more digits — left open.

## 4. Corner law status (task 3)

- d¹ law: **proved and universal** for this family (Section 2); the proof used only straightness of the p = 1 edges (S₀″ = −S₀) and corner/bisector boundary data, so it extends to general p-smoothed polygons Σℓ_j(x)^p = 1 — **conjecture:** P′(1) = Σ_corners (edge-length factor)·R₁(corner data); untested on non-regular examples (time ran out).
- d²/d³ corner law: **not found.** My attempt to split P″ into boundary + remainder integral by the same by-parts trick fails its numerical check (J₂+boundary − P″ = 5.26 at m=2 — the order-2 boundary bookkeeping is wrong as derived, and at order 3 the by-parts terms separately diverge like log θ and need regularisation). Reported honestly as open; the m=3 totals are nonetheless identified (Section 3).
- ζ(3) at every corner: m=2 PROVED (6√2), m=3 NUMERICAL (−9√3/8). Angle-dependence is nontrivial (sign change); m=4 value computed but unidentified. **CONJECTURE:** ζ(3) appears at every corner with an angle-dependent coefficient in Q(cos(π/2m)), with all other d³ constants elementary (logs of cyclotomic numbers, π², π²log, log³) — i.e. Clausen/log-sine values cancel at the diamond end, mirroring m=2.

## Evidence table

| Item | Check | Residual | Digits | Status |
|---|---|---|---|---|
| corner-arc rep. vs brute force (m=3,4; p≤3) | per-segment quad | 0 – 1.1e-50 | 50 | NUMERICAL |
| P(1) formula vs polygon perimeter (m=2..6) | closed form | ≤6.2e-61 | 60 | PROVED |
| m=2: P(1),P′,P″,P‴ vs R2 (PROVED there) | coefficient integrals | 6.2e-61, 2.1e-61, 2.5e-60, 5.5e-58 | 60,60,59,57 | PROVED (via R2) |
| m=3 P′, P‴ vs Richardson FD of direct P(p) | independent FD | ~5e-20, ~1e-11 | 19, 11 | NUMERICAL |
| P′(1) = P(1)R₁(0) (m=3,4) | coefficient integral | 4.2e-81 | 80 | PROVED (by-parts) + NUMERICAL |
| m=3 P′ = 2√3 log2 | PSLQ [1,−2] | 2.5e-90 | 90 | PROVED (corollary) |
| m=4 P′ = P(1)[log(1+√2)+(1−√2/2)log2] | PSLQ, maxcoeff 80 | 1.6e-90 | 90 | PROVED (corollary) |
| m=3 P″ closed form | PSLQ 25-basis, maxcoeff 32 | 9.9e-90 | 90 | NUMERICAL |
| m=3 P‴ closed form | PSLQ 32-basis, maxcoeff 96 | 1.9e-139 | 139 | NUMERICAL |
| m=4 P″ large-basis relation | maxcoeff 19354 | 2.4e-109 (of 170) | — | REJECTED (artifact risk) |
| J₂ by-parts decomposition | vs P″ totals | fails by O(1) | — | REJECTED/open |

## Scripts (all in `D15_polygon_p1/`; auto-saved)

`step_13.py` geometry + P(1) formula + first representation check; `step_14.py`–`step_15.py` brute-force per-segment + polygonal cross-checks; `step_16.py`–`step_20.py` coefficient-integral development (T₀′ fix, Q₀ sign fix — each bug caught by FD checks); `step_21.py`–`step_22.py` FD validation of ρ, Q expansions; `step_23.py` validated coefficients m=2..6; `step_24.py` 100-digit m=3,4 values (`vals_m34.json`); `step_25.py` P′(m=3) = 2√3log2 discovery; `step_26.py` P′(m=4) PSLQ; `step_27.py` P″(m=3) PSLQ; `step_28.py`/`step_29.py` 170-digit values (typo caught via FD; corrected run saved `vals_m34_170.json`); `step_30.py`–`step_31.py` P‴(m=3) PSLQ (Li₃(1/2) red herring removed); `step_32.py`–`step_33.py` m=4 P″ attempts (rejected); final J₂ check (failed, reported open).

**Bottom line:** at the diamond end, ζ(3) does appear at the hexagon corners (m=3) with explicit coefficient −9√3/8 (and at m=2 with 6√2, proved in R2); P′(1) obeys the clean universal law P(1)·R₁(0) (proved), while P″ loses the pure-π² form (logs enter) and P‴ stays in the elementary+ζ(3) universe with Clausen/log-sine constants cancelling — m=4 identification and the general d²/d³ corner law remain open, with 170-digit targets ready in `vals_m34_170.json`.