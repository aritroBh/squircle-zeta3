# AV5c — Square-type ends rechecked with the corner-polylogarithm bases

**Verdict.** AV5's "no zeta at orders 1–4" at the squircle's square end is **wrong**, and for the reason diagnosed: its PSLQ basis (step_57: `{1, π, π², π³, ln2, ln²2, ln³2, Catalan, ζ(3), πln2, π²ln2, π·C, √2, √2π, √2ln2}`) contained **no** `L = ln(1+√2)`, no `Li₂(b)`, no `Li₃(±b)`, no `Li₃(1/√2)` — i.e. none of the corner alphabet. With the correct basis the same coefficients identify immediately: **weight 1 (logs) at order ε¹, π² + Li₂(√2−1) at ε², ζ(3) + L(χ₈,3) + Li₃(b²) + Li₃(1/√2) at ε³, tetralogarithms at ε⁴** — all PROVED by the paper's thm:M1, thm:M2, thm:W4R, independently recomputed here, and rediscovered blind by LLL at two precisions. The same refit was carried out for the octagon end (done), partially for the hexagon and Lamé ends (status below).

Notation: `u = 1/p`, `ε = q−1 = 1/(p−1) = u/(1−u)` (AV5's squircle parameter; note AV5's hexagon/octagon/Lamé fits used `ε = 1/p = u` directly), `a = 1+√2`, `b = √2−1 = 1/a`, `L = ln a`, `λ = ln2`, `c = cos(α/2)` for corner angle α, `χ₃(x) = (Li₃(x)−Li₃(−x))/2`.

---

## 1. Control: the squircle's square end |x|^p + |y|^p = 1, p → ∞

**Convention dictionary (exact).** Paper: `P(p) = 8(1 + M₀/p + M₁/p² + M₂/p³ + M₃/p⁴ + …)`; AV5: `P = 8 + d₁ε + d₂ε² + d₃ε³ + d₄ε⁴ + …`. Since `ε = u/(1−u)`,
`d₁ = 8M₀`, `d₂ = 8(M₁−M₀)`, `d₃ = 8(M₂−2M₁+M₀)`, `d₄ = 8(M₃−3M₂+3M₁−M₀)`.
All four reproduce AV5's fitted numbers to AV5's own fit precision (27, 19, 13, 13 digits respectively) — AV5's *numbers* were right; only its identification failed.

**Identifications** (λ = ln2, L = ln(1+√2), b = √2−1):

| order | coefficient | identified constants | status |
|---|---|---|---|
| ε¹ | d₁ = −4.42642639776412573896646431234367549780712… | `8λ − 8√2 L` (weight 1, logs only) | **PROVED** (Ahmed; paper §M₀) |
| ε² | d₂ = 6.44737916742048409591337248006891438207614… | `(19√2/6−1/3)π² + 4λ² + 16√2λL − (8+20√2)L² − 40√2 Li₂(b) − 8λ + 8√2 L` — contains the dilogarithm **Li₂(√2−1)** (weight 2 beyond products) | **PROVED** (thm:M1) + LLL-rediscovered |
| ε³ | d₃ = −8.64562722577499764825032014971248538616879… | **ζ(3) with coefficient `189√2/8 − 19/2 = 23.9107954110643705279398961…`**, plus **320·L(χ₈,3) + 22·Li₃(b²) − 216√2·Li₃(1/√2)** plus products (full formula below) | **PROVED** (thm:M2) + LLL-rediscovered |
| ε⁴ | d₄ = 11.20333474476195901563115139921301091208657… | Li₄ at `b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2` and ζ(3)×log products (via M₃ = R+W₃′+T) | **PROVED** (thm:W4R) + NUMERICAL (residual 1.3e-99) |

The d₃ formula (verified to 2.0e-238, raw Li₃(±b) form to 4.2e-240, both at 250 dps):
```
d₃ = (189√2/8 − 19/2)ζ(3) + 320·L(χ₈,3) + 22·Li₃(b²) − 216√2·Li₃(1/√2)
     + (80√2 + 96L − 28√2λ)Li₂(b) + (1 − 19√2/3)π² − (14/3)π²L − (1/3 + 65√2/6)π²λ
     − 6λ² + (4/3 + 9√2/2)λ³ + 11√2λ²L − 32√2λL − (8+14√2)λL²
     + (16+40√2)L² + (100/3 + 4√2/3)L³ + 8λ − 8√2L ,
```
where `22 Li₃(b²) + 320 L(χ₈,3)` is `(88+120√2)Li₃(b) + (88−120√2)Li₃(−b)` after duplication and the paper's proved bridge `48χ₃(b) = 32√2 L(χ₈,3) − 3π²L + 4L³` (bridge re-verified here to 2.8e-250).

**Evidence (all printed):**
- *Direct perimeter fit* (independent of all formulas): P(p) at 160 dps, p = 32…96, 6-term solve, validation residuals ≤ 3.8e-18: recovers M₀ to 13 digits, M₁ 11, M₂ 9, M₃ 7 (fit-noise limited; earlier p = 16…48 run gave 11/9/7/6). NUMERICAL.
- *Corner-integral quadrature* (independent): `M(u;π/2) = ∫₀¹[(1+g)^{-1-u}√(1+g^{2−2u}) − 1]g^{u−1}dg` expanded in u with sympy and integrated term-by-term (110 dps): M₀…M₃ match the closed forms to ~100 digits (residuals 2.3e-101, 4.1e-100, 8.4e-100, 1.3e-99). NUMERICAL.
- *d₂ LLL/PSLQ rediscovery*, basis `{π², λ², L², λL, Li₂(b), λ, L, 1}×{1,√2}`: **identical relation at dps 160 and dps 210, residual 0.0 at 235 dps, maxc 240** — exactly the thm:M1 content. NUMERICAL(2 precisions) + PROVED by thm:M1.
- *d₃ LLL rediscovery*, basis `{ζ(3), L(χ₈,3), Li₃(b²), Li₃(1/√2), π²λ, π²L, λ³, λ²L, λL², L³, λLi₂(b), LLi₂(b), π², λ², L², λL, Li₂(b), λ, L, 1}×{1,√2}` (40 elements): **identical 27-term relation at C = 1e195/dps 215 and C = 1e165/dps 185, maxc 7680, residual 5.8e-233 at 235 dps** — matches the formula above coefficient-by-coefficient. NUMERICAL(2 precisions) + PROVED by thm:M2.
- *d₃ without trilogarithms* (basis `{ζ(3), products, w₂, w₁}×{1,√2}`, 34 elements): only spurious candidates at C = 1e180 and C = 1e230; their residuals at 245 dps are 4.0e-176 and 5.1e-196 (not exact) with coefficients growing 7e4 → 5e5 as C grows — **no exact relation**: trilogarithms are genuinely needed (NUMERICAL, mirroring the paper's 290-digit statement for M₂).
- *AV5's own basis*: d₂ and d₃ give **NONE** with maxcoeff ≤ 10⁶ at C = 1e150/dps 170 and C = 1e210/dps 230 — this reproduces AV5's false negative and identifies its cause.

So at the square end, **ζ(3) first appears at cubic order** (ε³, equivalently u³ — same coefficient `189√2/8 − 19/2` in both conventions, since only M₂ carries ζ(3)), together with three more weight-3 constants; the weight rises by one each order: `wt(d_k) = k`.

---

## 2. The corner dictionary used for the other ends

Near a polygon corner where exactly two edge-functionals are active, the p-norm curve is the paper's oblique elbow `x^p + y^p = 1` mapped linearly by the two edge vectors `e₁, e₂` (inactive terms are exponentially small); the perimeter deficit of the corner is
`δ(e₁,e₂;u) = ∫_{fan}[((c·e₁)^q + (c·e₂)^q)^{1/q} − c·(e₁+e₂)]dθ`, `q = 1/(1−u)`, fan = `{c·e₁, c·e₂ ≥ 0}`.
For unit-speed edges the paper gives `δ = 2u·M(u;α)`, `M(u;α) = Σ M_k(α)u^k`, and my term-by-term quadrature of this integral (sympy u-expansion + 110/300-dps tanh-sinh) reproduces: the elementary `M₀(α) = −(1+c)ln(1+c) − (1−c)ln(1−c)` (3 angles, residuals ~1e-101 and 3.4e-281), the paper's proved `M₁(2π/3)` formula (residual 1.8e-101), and all of the squircle's M₀…M₃ (§1).

---

## 3. Octagon end |x|^p + |y|^p + 2^{−p/2}(|x+y|^p + |x−y|^p) = 1 (eight 135° corners, edge speeds (√2,√2))

Prediction: `P_oct(p) = 16(√2−1) + 16√2·Σ_{k≥0} M_k(3π/4)·u^{k+1}`.

**Validation against AV5's direct perimeters** (their ε = u): differences `7.4e-11, 7.1e-12, 6.8e-13, 6.4e-14, 6.1e-15, 5.8e-16, 5.6e-17` at p = 50…839, decaying like the omitted u⁵ term — the corner decomposition is NUMERICALLY exact to AV5's data; predicted coefficients match AV5's fits to 19/12/8 digits (AV5's fit precision).

| order (u) | coefficient | constants | status |
|---|---|---|---|
| u¹ | `16√2 M₀(3π/4)` = −3.39973594013030381743550157… | `16√2·[−(1+c)ln(1+c) − (1−c)ln(1−c)]`, c = cos 3π/8: logs only (in the degree-4 field Q(c) ∋ ln(1±c)) | **PROVED** (elementary; quad 3.4e-281) |
| u² | `16√2 M₁(3π/4)` = 1.66330853474362645505407835… | dilogarithms (weight 2): proved 27-dilog form (thm:Malpha; matches 300-dps quad to 1.3e-282) and the LLL short form below | **PROVED** (formula) / NUMERICAL (short form) |
| u³ | `16√2 M₂(3π/4)` = −0.12157940917958534249642464… | weight-3 **trilogarithms** (thm:M2alpha: `ζ(2)Φ(c) − ½Li₂(c²) + Y(c) + Y(−c)`, HPLs at ±c; **no ζ(3) as a function of the angle**); 280-digit value computed; short-form LLL not finished | value NUMERICAL (280 digits); structure **PROVED†** (paper); closed form OPEN |
| u⁴ | `16√2 M₃(3π/4)` = 0.13780491163662157980552077… | **−ζ(3)·Φ(c)** with `Φ(c) = ½ − (1−c²)/(2c)·artanh c ≈ 0.0509 ≠ 0`, plus `ζ(2)Ψ(c)` and a Q[[c²]]-series (thm:flat) | **PROVED†** (paper, modulo standard 1,π²,ζ(3)-independence caveat); my independent coefficient check OPEN |

M₁(3π/4) short form (LLL, 51-dim lattice, C = 1e260/dps 280, norm² = 369, maxc 11, **residual 2.7e-282 at 290 dps**; consistent with the paper's "reduced form checked to 500 digits"):
```
2M₁(3π/4) = −cλ² + 5cλL − 2ln(1−c)² − 4c·ln(1−c)² − 3λln(1−c) − 11cλln(1−c) + 2Lln(1−c)
            + 10cLln(1−c) − 4c·Li₂(r) + 5c·Li₂(b²) − Li₂(c²) + 4c·Li₂(c²) − 4c·Li₂(c)
            + 3c·Li₂(1/2) + c·Li₂(1−1/√2) + c·Li₂(2−√2) ,   r = (1−c)/(1+c) = tan²(3π/16).
```
En route the LLL produced the exact dilog identities `Li₂(c²) = λln(1−c²) − Li₂(b²) − Li₂(1/2) + Li₂(1−1/√2) + Li₂(2−√2)`, `Li₂(s²) = Li₂(c²) + 3Li₂(b²) + L² + 2L·ln(s²)` (s = sin 3π/8), residuals ≤ 2e-236 — the corner alphabet at 135° really is the `c = cos 3π/8` field, as the paper's field rule says.

**Octagon verdict:** AV5's "period constants, no zeta" is wrong in kind: dilogarithms at u² (PROVED), trilogarithms (no ζ(3)) at u³ (PROVED† structure; value computed), **ζ(3) at u⁴** with exact coefficient `−16√2·Φ(cos 3π/8)` (PROVED†).

---

## 4. Hexagon end |x|^p + |y|^p + |x+y|^p = 1

Geometry: corners at (±1,0),(0,±1) are **four 135° corners with mixed edge speeds (√2, 1)**; corners at (±1,∓1) are **two 90° corners with unit speeds — literally the squircle's elbow**. Hence
`P_hex(p) = (4 + 2√2) + 4·δ(√2,1,3π/4;u) + 4u·M(u;π/2) + O(e^{−cp})`,
`δ = ∫₀^{π/4}[((cosθ−sinθ)^q + sin^qθ)^{1/q} − cosθ]dθ` (derived here; its leading term δ₁ is elementary: the log-trig integrals over [0,π/4] reduce to `λ, L, √2` — no Catalan).

- The two right corners alone contribute `4M₂(π/2)·u³` at cubic order, i.e. the squircle's **ζ(3) + L(χ₈,3) + Li₃(b²) + Li₃(1/√2)** combination with coefficient `4·(189√2/64 − 19/16) = 189√2/16 − 19/4` on ζ(3): **weight-3 constants provably enter the hexagon end at u³** (PROVED, via the corner decomposition + thm:M2).
- The total u²/u³ coefficients = `4δ_k + 4M_{k−1}(π/2)`: the mixed-corner δ₂, δ₃ (their own 135°-field polylogs, possibly trading ζ(3)) were not finished computing — the *total* identification is **OPEN**. AV5's fitted d₁ = −3.635349946172…, d₂ = 1.718558284652…, d₃ = −0.138363729458… (their u-convention) remain NUMERICAL; their "no zeta" claim is superseded by the proved ζ(3) content above.

## 5. Lamé ends |x|^p + |y|^{λp} = 1, p → ∞ (λ = 2, 3, 1/2)

- The support function satisfies the exact implicit system `A^{1/(1−u)} + B^{λ/(λ−u)} = 1`, `h = c₁A^{u/(1−u)} + c₂B^{u/(λ−u)}`, `A = c₁/z`, `B = c₂/(λz)` (derived here). Its first order gives AV5's profile `F(θ) = c₁ln(λc₁) + (c₂/λ)ln c₂ − (c₁+c₂/λ)ln(λc₁+c₂)` exactly, hence d₁ = 4∫F dθ in AV5's closed form — **PROVED** (elementary: `ln λ, ln2, ln(1+λ²), ln tan(π/8 + arctan λ /2)`-type; at λ = 1 it reduces exactly to the squircle's `8λ − 8√2 L`, checked symbolically).
- d₂, d₃ (AV5's numbers, ~15 reliable digits, u-convention): **OPEN** — the higher-order recomputation and the LLL against the λ-adapted basis (natural argument `t_λ = λ/(1+√(1+λ²)) = tan(arctan λ /2)`: `t₂ = (√5−1)/2`, `t₃ = 3/(1+√10)`; note `t₁ = b`) were not completed. AV5's "no zeta" here is unreliable (its basis had no polylogs) and must be redone.

## 6. AV5 negatives that survive

Soft-max square and LSE diamond: `P = P₀ + c₁/p + O(e^{−2p})` with elementary c₁ — there are *no* algebraic orders beyond 1/p, so nothing can hide; Minkowski sum (`8 + 2πr` exact) and the sech family (`8 + 4πw² + O(e^{−1/w})`) similarly. These four negatives stand (PROVED/NUMERICAL as in AV5). The Lamé p → 1⁺ end (exponentially narrow physical layer, analytic bulk) is not a square-type end and was not rechecked.

## 7. Corrected rules

**R1 (unchanged, correct).** A *shrinking* softplus layer (two flat fans of the support function meeting with algebraic overshoot, exponential tails: profile `log(1+e^{−a|s|})` embedded as `h = m·e^{τG}`) gives η-moments: ζ(n+1) at order τ^{n+1}, single-zeta collapse through weight 7, alternating MZVs at weight 8. (Squircle diamond end; untouched by this recheck.)

**R2′ (corrected R2).** A *non-shrinking* corner profile does **not** mean "period constants only". At any polygon end smoothed by the p-norm, each corner contributes the corner series `M(u;α)` (unit-speed) or its mixed-speed analogue, and the coefficient at order u^k (equivalently the leading behavior in ε = 1/(p−1)) has **polylogarithmic weight k**: logs at k = 1, π² + dilogs at k = 2, **trilogarithms and ζ(3) at k = 3**, Li₄-class constants at k = 4. The alphabet is the corner field's: `b = √2−1, 1/√2, 1/2` at α = π/2; S-units of `Q(cos(α/2), sin(α/2))` in general (paper's field rule, thm:conj:field). ζ(3) itself appears provably from weight 3 on: at a right angle already in M₂ (thm:M2); at a general angle first in M₃(α) with the explicit coefficient `−Φ(cos α/2)` (thm:flat). What survives of the old R2: the *mechanism* is different from R1 — no softplus moments, no η-values; the zeta values enter through corner iterated integrals, not through a shrinking layer.

**R3 (unchanged).** A symmetric/rigid/even profile (log cosh, sech, flat shift) has vanishing odd moments: no ζ(3) (soft-max square, LSE diamond, Minkowski, sech family — all confirmed).

## 8. Plain words

AV5 asked the right question — "is the squircle's ζ(3) universal?" — but looked for the answer with the wrong dictionary. At the square end of the squircle the perimeter coefficients are not made of π, ln 2 and Catalan at all: they are built from logarithms of 1+√2 and from polylogarithms at √2−1 and its relatives. Since none of those constants were in AV5's search basis, the lattice reduction returned "no relation", which AV5 read as "no zeta". In fact ζ(3) was sitting in the very next coefficient it printed: d₃ contains ζ(3) with a coefficient near 24, together with three independent trilogarithmic constants, and this is not numerology — it is the paper's Theorem M₂, proved there and re-verified here three independent ways (direct high-precision perimeter fit, term-by-term corner quadrature, and blind LLL rediscovery at two precisions). The corrected picture is nicer than the old one: *every* p-norm polygon end carries the corner's polylogarithm ladder, one weight per order; the angle of the corner only changes which field the constants live in. The square end is not an exception to the squircle's ζ(3) story — it is the same story told in a different accent.

## 9. Evidence summary (residuals and digits, all printed)

| claim | evidence | digits/residual | status |
|---|---|---|---|
| d₁ = 8λ − 8√2L | closed form vs AV5 fit | 27 digits | PROVED |
| d₂ formula (thm:M1) | exact relation at 250 dps | 8.4e-241 | PROVED |
| d₂ LLL rediscovery | identical at dps 160 & 210, 16-elt basis | 0.0 @235 dps, maxc 240 | NUMERICAL×2 (+PROVED) |
| d₃ formula (thm:M2) | raw & reduced forms at 250 dps | 4.2e-240 / 2.0e-238 | PROVED |
| d₃ LLL rediscovery | identical at C=1e195/1e165, 40-elt basis | 5.8e-233 @235 dps, maxc 7680 | NUMERICAL×2 (+PROVED) |
| d₃ needs trilogs | ζ(3)+products basis, C=1e180/1e230 | spurious only, 4.0e-176/5.1e-196 | NUMERICAL |
| d₂,d₃ in AV5 basis | C=1e150/1e210, maxcoeff 1e6 | none | NUMERICAL (explains false negative) |
| d₄ = 8(M₃−3M₂+3M₁−M₀) | thm:W4R closed form; corner quad | 1.3e-99 | PROVED + NUMERICAL |
| direct P(p) fit (control) | p=32…96, 160 dps | M₀…M₃ to 13/11/9/7 digits | NUMERICAL |
| χ₃(b) bridge | eq. (chi3) at 250 dps | 2.8e-250 | PROVED (paper) + check |
| M₀(α) elementary, 3 angles | quad vs formula | ~1e-101; 3.4e-281 (3π/4) | PROVED + NUMERICAL |
| M₁(2π/3) paper formula | quad vs formula | 1.8e-101 | PROVED† + NUMERICAL |
| octagon = 16√2·uM(u;3π/4) | vs AV5 direct P(p), p=50…839 | 7.4e-11→5.6e-17 (u⁵ tail) | NUMERICAL (decomposition) |
| M₁(3π/4) 27-dilog formula | thm:Malpha vs 300-dps quad | 1.3e-282 | PROVED + NUMERICAL |
| M₁(3π/4) short form | LLL C=1e260/dps280 | 2.7e-282, maxc 11 | NUMERICAL |
| octagon ζ(3) at u⁴ | `−16√2·Φ(cos 3π/8)·ζ(3)` | Φ ≈ 0.0509 ≠ 0 | PROVED† (paper) |
| hexagon ζ(3) at u³ | 2 square corners contribute `4M₂(π/2)` | — | PROVED (component); total OPEN |
| Lamé d₁ (λ=2,3,1/2) | profile re-derived = AV5's F; λ=1 ⇒ squircle d₁ | exact | PROVED |
| Lamé d₂,d₃; hexagon totals; M₂(3π/4), M₃(3π/4) closed forms | not completed | — | OPEN |

## 10. Files (work dir `new_results/AV5c_square_ends_recheck`)

- `progress.log` — one line per finished step (steps 29–64).
- `step_29.py` constants module (77 constants, 250 dps; χ₃ bridge check) · `step_30.py` paper closed forms M₀–M₃ & binomial d₁–d₄ · `step_31.py` sympy u-expansion of the corner integrand · `step_32.py` M_k(α) quadrature (π/2, 3π/4, 2π/3) · `step_33.py` cross-checks + d₁–d₄ vs closed forms · `step_34–36.py` independent direct perimeter fits · `step_37–38.py` exact d₂/d₃ relations at 250 dps · `step_39–42.py` d₂ PSLQ (spurious-control + clean 160/210-dps rediscovery) · `step_43–48.py` d₃ LLL rediscovery (two precisions) · `step_50–52.py` AV5-basis and no-trilog controls · `step_53–54.py` M₁(3π/4) formula check + octagon validation vs AV5 data · `step_55–59.py` 135°-basis reduction (exact dilog identities) · `step_60–61.py` M_k(3π/4) at 300 dps · `step_62–64.py` M₁(3π/4) 282-digit check and LLL short form. (step_01–28: reading/derivation phase; step_44 was the flint-LLL prototype.)
- Checkpoints: `ckpt_constants.json`, `ckpt_paper_Mk.json` (M₀–M₃, d₁–d₄ at 240 digits), `ckpt_corner_integrand.pkl`, `ckpt_Mk_cornerintegral.json`, `ckpt_sq_direct_Ps*.json`, `ckpt_pslq_d2.json`, `ckpt_pslq_d3.json`, `ckpt_Mk_3pi4_hi.json` (280-digit M_k(3π/4)), `ckpt_M1_3pi4_formula.json`.

**Bottom line for the corrected AV5 table:**

| end | order | constants found | digits | basis | status |
|---|---|---|---|---|---|
| squircle p→∞ | ε¹ | λ, ln(1+√2) | 250 | — | PROVED |
| squircle p→∞ | ε² | + π², **Li₂(√2−1)** | 250; LLL×2 (0.0 @235) | {π²,λ²,L²,λL,Li₂(b),λ,L,1}×{1,√2} | PROVED + NUMERICAL |
| squircle p→∞ | ε³ | + **ζ(3) (coeff 189√2/8−19/2), L(χ₈,3), Li₃(b²), Li₃(1/√2)** | 250; LLL×2 (5.8e-233) | +{ζ(3),Lχ,Li₃(b²),Li₃(1/√2),w3-products} | PROVED + NUMERICAL |
| squircle p→∞ | ε⁴ | + Li₄(b,−b,−b²,1/√2,1/2,1−1/√2,2−√2) | 250; quad 1.3e-99 | — | PROVED + NUMERICAL |
| octagon p→∞ | u¹ | logs of 1±cos(3π/8) | 280 | — | PROVED |
| octagon p→∞ | u² | + π² and **dilogs (Li₂(c²), Li₂(r), Li₂(b²), …)**, c=cos 3π/8 | 282; LLL 2.7e-282 | 135° basis ×{1,c} | PROVED (formula) + NUMERICAL |
| octagon p→∞ | u³ | **trilogarithms (HPLs at ±c); no ζ(3)** | 280 | thm:M2alpha | PROVED† + NUMERICAL value |
| octagon p→∞ | u⁴ | **ζ(3)** with coeff −16√2·Φ(cos 3π/8) | — | thm:flat | PROVED† |
| hexagon p→∞ | u³ | **ζ(3), L(χ₈,3), Li₃(b²), Li₃(1/√2)** from its two square corners (coeff on ζ(3): 189√2/16−19/4) + 135°-corner content | — | — | PROVED (component); total OPEN |
| Lamé p→∞ (λ=2,3,1/2) | u¹ | elementary (ln λ, ln(1+λ²), ln tan) | exact | — | PROVED |
| Lamé p→∞ | u², u³ | unknown (AV5's "no ζ" unreliable — no polylogs in basis) | — | t_λ-basis | OPEN |

Corrected rules: **R1** (shrinking softplus layer ⇒ ζ(n+1) at τ^{n+1}) and **R3** (even/rigid profile ⇒ no ζ(3)) stand; **R2′**: non-shrinking p-norm corners bring in the corner-angle polylog ladder — weight k at order k, with ζ(3) and trilogarithms arriving at cubic order (right angle: thm:M2; general angle: trilogs at u³, ζ(3) explicitly at u⁴ via −Φ(cos α/2)).