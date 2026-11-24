# AV5: How universal is the squircle's ζ(3) mechanism?

**Headline.** The squircle's ζ(3)-at-τ³ mechanism is **not universal** — it is exceptional. Across ten other polygon-degenerating families/ends, ζ(3) never appears. The deciding feature is the *layer profile of the support-function deviation at the polygon's edge normals*: only when two flat fans of the polygon's support function meet and the curve overshoots the tent algebraically does the **log-sum-exp (softplus) profile** log(1+e^{−2|s|}) form, whose moments are η-values (1−2^{1−n})ζ(n) — so ζ(2), ζ(3), ζ(4), … appear at orders τ², τ³, τ⁴, …. All other ends have either a *non-shrinking* bump (period constants, no zeta), an exponentially narrow layer (analytic constants, no zeta), or a symmetric/rigid profile (no odd moments, no ζ(3)).

---

## 0. Calibration: the squircle |x|^p+|y|^p = 1 (mechanism confirmed and completed)

### Diamond end (p → 1⁺, τ = 1 − 1/p)
**Perimeter expansion** (P = 8∫₀^{π/4}(c^q+s^q)^{1/q}dθ, q = 1/τ; computed to 75–90 digits and fitted by geometric-τ interpolation):

| coeff | value | identification | status |
|---|---|---|---|
| c₂ | 2.32628806654629317805252031205646842997329 | √2 ζ(2) | **PROVED** |
| c₃ | 1.69996517519250036777674394595165573772495 | √2 ζ(3) | **PROVED** |
| c₄ | −3.15693715422070561248491036997352611509387 | −(33/16)√2 ζ(4) | NUMERICAL* |
| c₅ | −5.99340657187070660107725941712303200938264 | −(√2/4)(7ζ(2)ζ(3)+3ζ(5)) | NUMERICAL* |
| c₆ | 20.9114781464325920489424555553583910786688 | (√2/160)(577ζ(2)³−140ζ(3)²) | NUMERICAL* |
| c₇ | 63.544518077654247431367031506734473687014 | (√2/160)(1071ζ(2)²ζ(3)+1410ζ(2)ζ(5)+1290ζ(7)) | NUMERICAL* |
| c₈ | −300.899292485393445660523299799385402718211 | not in single-zeta algebra | NUMERICAL (OPEN formula) |
| c₉ | −1287.73316060792062005836415215762506093905 | — | NUMERICAL (OPEN formula) |

\* assembly exact; constants identified by PSLQ at two precisions (dps 60 & 120), residuals < 10⁻¹¹⁰, pigeonhole margins 10⁸–10¹⁰⁰.

**Layer profile** (derived symbolically with sympy, verified numerically): with m(θ) = max(cos θ, sin θ) the diamond's tent,
h(π/4+τs) = m·exp(τ g₁ + τ³ g₃ + τ⁵ g₅ + τ⁷ g₇ + …), where
g₁(s) = log(1+e^{−2|s|}) (even — **the softplus profile**), g₃ = (2/3)|s|³(tanh|s|−1), g₅ = (2/3)|s|⁵(tanh|s|−1) + (2/9)s⁶sech²s, g₇ = −(8/81)|s|⁹tanh s sech²s + (4/9)s⁸sech²s + (244/315)|s|⁷(tanh|s|−1). Even coefficients g₂=g₄=g₆=0.

Moments: ∫₀^∞ sʲ g₁ ds = j! η(j+2)/2^{j+1} (j=0: ζ(2)/4; j=1: 3ζ(3)/16; j=2: 7ζ(4)/32); ∫₀^∞ g₁² ds = ζ(3)/8; ∫ sᵃ g₁ᵇ ds = colored MZV of weight a+b+1. Coefficient c_{k+1} = 4√2∫₀^∞ f_k(s) ds with f_k universal polynomials in (g₁,g₃,g₅,g₇,s). ζ(3) at τ³ comes from ∫(g₁²/2 + s g₁)ds = ζ(3)/16 + 3ζ(3)/16 — i.e., **the square of the profile plus the tent's |s|-slope times the profile** (the asymmetric nesting).

**The single-zeta collapse**: the assembled coefficients c₂…c₇ (weights 2–7) lie entirely in the Riemann-zeta algebra — every ln 2, Li_n(1/2), and alternating constant in the profile moments *cancels* (F₃ = −(33/160)ζ(2)², F₄ = −(7ζ(2)ζ(3)+3ζ(5))/16, F₅ = (577ζ(2)³−140ζ(3)²)/640, F₆ = (1071ζ(2)²ζ(3)+1410ζ(2)ζ(5)+1290ζ(7))/640; residuals 6.2e-114, 3.3e-114, 2.9e-112, 3.4e-113 at dps 125). The collapse **breaks at weight 8**: c₈ is not a zeta monomial (PSLQ at dps 200, maxcoeff 10¹²).

**Fit-vs-prediction** (direct perimeter vs layer-moment coefficients): 26.1, 21.0, 16.7, 12.9, 9.5, 6.5, 3.9, 1.7 agreeing digits for c₂…c₉ — the layer expansion is the true asymptotic series. Sign pattern ++−−++−− (period 4).

**Area (PROVED, exact)**: A(τ) = 2·exp(Σ a_k τ^k), a_k = (1/k)[1 + (2−2^k)(ζ(k)−1)]: a₁ = 1, a₂ = 3/2−ζ(2), a₃ = 7/3−2ζ(3), a₄ = 15/4−(7/2)ζ(4), a₅ = 31/5−6ζ(5). ζ(3) enters at τ³ here too (coefficient −2 in the log).

### Square end (p → ∞, ε = q−1 = 1/(p−1))
P = 8 + d₁ε + d₂ε² + d₃ε³ + d₄ε⁴ + … with d₁ = −4.426426397764125738966464267, d₂ = 6.447379167420484096, d₃ = −8.645627225775, d₄ = 11.203334744762. **Non-shrinking profile** g(θ) = cos θ[u ln u − (1+u)ln(1+u)], u = tan θ; d₁ = 8∫₀^{π/4}g dθ verified to 28 digits. Constants: NUMERICAL, no zeta at orders 1–4 (PSLQ none, two precisions). Layer at the vanishing branch is exponentially narrow (~e^{−1/ε}).

---

## 1. Anisotropic Lamé curves |x|^p + |y|^{λp} = 1 (λ = 2, 3, 1/2)

Areas are exact: A = 4Γ(1+1/p)Γ(1+1/(λp))/Γ(1+1/p+1/(λp)).

### p → 1⁺ (limit: digon |x|+|y|^λ = 1, corners at (0,±1))
P(1,2) = P(1,1/2) = 2√5 + asinh(2) = 2√5 + ln(2+√5) (**PROVED exact**); P(1,3) = 6.191462618734441 (NUMERICAL).

| λ | c₁ (τ) | c₂ (τ²) | c₃ (τ³) | status |
|---|---|---|---|---|
| 2 | 1.003341875013133417 | 0.281435993067217 | −1.752129907993865 | NUMERICAL, no ζ |
| 3 | 1.139693253046047103 | −0.338361698717748 | −0.680245366482988 | NUMERICAL, no ζ |
| 1/2 | −1.222923907410829532 | 2.163782762941897 | −0.563922871331084 | NUMERICAL, no ζ |

Physical corner layer is **exponentially narrow**: δh inside the corner fan decays super-algebraically (2.7e-5 → 5.0e-8 → 3.1e-13 → 2.4e-28 as τ halves), while the bulk deforms analytically at O(τ) (δh/τ → 0.408). Hence an analytic expansion with period constants — PSLQ finds no ζ(2)/ζ(3) at any coefficient (two precisions, bases up to 14 elements).

### p → ∞ (limit: square)
**Non-shrinking profile** (derived exactly): F(θ) = cos θ ln(λ cos θ) + (sin θ/λ) ln sin θ − (cos θ + sin θ/λ) ln(sin θ + λ cos θ), and P = 8 + d₁ε + … with **d₁ = 4∫₀^{π/2}F dθ in closed form (PROVED)**:
d₁ = 4[ln λ + (ln 2 − 1)(1 + 1/λ) − I(λ)/λ], I(λ) = √(1+λ²)[(sinφ+cosφ)ln√(1+λ²) + G(φ+π/2)−G(φ)], φ = arctan λ, G(ψ) = −cosψ ln sinψ + ln tan(ψ/2) + cosψ.

| λ | d₁ (closed form) | d₂ | d₃ | status |
|---|---|---|---|---|
| 2 | −3.063001374800518 | 0.969769008528502 | −0.057175053431166 | d₁ PROVED; rest NUMERICAL, no ζ |
| 3 | −2.421476812773221 | 0.607958095070304 | −0.026980951795464 | same |
| 1/2 | −6.126002749601054 | 3.879076034135909 | −0.457400437755694 | d₁ = 2d₁(2) (exact λ↔1/λ symmetry) |

Closed form vs fit: 11 digits agreement (fit precision limit). Elementary constants only (λ, ln 2, ln(1+λ²), ln tan) — no zeta.

---

## 2. Soft-max square and log-sum-exp diamond

**Soft-max square** log(cosh px + cosh py)/p = 1. Support function **exact to all algebraic orders** (error O(e^{−2p}/p)): on |β| ≤ π/4,
h(π/4+β) = √2 cos β + (√2/p)F(β), F(β) = ½[(c+s)ln(c+s) + (c−s)ln(c−s)] − c ln c (even; verified vs direct maximization, ≤ 5e-10 at p=10).
Perimeter: **P = 8 + c₁/p + O(e^{−2p})** with c₁ = 8√2∫₀^{π/4}F dβ = **16 ln 2 − 8√2 ln(1+√2) = 1.118751046715436736371392659321737046796877152219** (PROVED elementary). No 1/p², 1/p³ terms (fit residuals at the e^{−2p} level), no zeta anywhere. Physical corner profile σ(δ) = 1 − (1/p) log cosh(pδ) — the even "1/cosh"-type profile.

**LSE diamond** (same family, rotated): P = 4√2 + (4√2 ln 2 − 8 ln(1+√2))/p + O(e^{−2p}) — verified to −2.8e-17 (p=20), −2.1e-28 (p=40). No zeta.

---

## 3. p-norm polygons: hexagon and octagon

**Hexagon** |x|^p+|y|^p+|x+y|^p = 1 (limit perimeter 2√2+4 = 6.828427; support function via Legendre dual h(c) = min_{z₃}(|c₁−z₃|^q+|c₂−z₃|^q+|z₃|^q)^{1/q}, validated at p=2 against the ellipse to 65 digits):
P = P_hex − d₁ε − d₂ε² − …, d₁ = −3.635349946172092407738869, d₂ = 1.718558284652393, d₃ = −0.138363729458309, d₄ = 0.148056777504145. Leading term **linear** in ε: edge normals pinned (h(π/4) = hexagon value to e^{−cp}), deviation a non-shrinking negative even bump Φ(θ) (zero at corner normals; Φ = (h_p−h_Π)/(q−1) converges: −0.4855, −0.6253, −0.5271 at θ = 0.196, 0.393, 0.589). NUMERICAL, no zeta (PSLQ, two precisions).

**Octagon** |x|^p+|y|^p+2^{−p/2}(|x+y|^p+|x−y|^p) = 1 (exact 45°-rotation symmetry — the norm itself is octagon-symmetric; limit 16(√2−1) = 6.627417): d₁ = −3.399735940130302868833826, d₂ = 1.663308534742409, d₃ = −0.121579408610481. Same structure: pinned, non-shrinking, NUMERICAL, no zeta.

---

## 4. Symmetric-layer families

**(a) Minkowski sum** (square ⊕ disk): h = h_sq + r, P = 8 + 2πr exactly (Steiner). No layer, no zeta. ✓ parity
**(b) sech-profile family** h_w = h_sq + w Σⱼ sech((θ−θ_j)/w) (valid convex support function: sech + sech″ = 2 sech tanh² ≥ 0): **P = 8 + 4π w² + O(e^{−1/w})** — verified to 0.0–1.4e-20 across w ∈ [0.0026, 0.02]. Single w² term (∫ sech = π); no w³, no ζ(3). ✓ parity confirmed: a rigid even profile contributes only its zeroth moment.

---

## 5. The rule, and its tests

**R1 (softplus-layer end → zeta values).** If at an edge normal θ_e of the polygon the support function overshoots the tent algebraically (h_τ(θ_e) − h_Π(θ_e) = cτ, c > 0 — bulging edges) and the deviation decays exponentially outside |θ−θ_e| ≲ τ, then the layer profile is φ(s) = log(1+e^{−a|s|}), embedded exponentially (h = m·e^{τG}); its moments are j!(1−2^{1−(j+2)})ζ(j+2) and colored MZVs; the coefficient of τ^{n+1} carries ζ(n+1): ζ(2) at τ², **ζ(3) at τ³** (from ∫φ²/2 and ∫|s|φ), ζ(4) at τ⁴, …; coefficients stay in the single-zeta algebra through weight 7, alternating irreducibles enter at weight 8.
**R2 (non-shrinking or pinned profile → period constants only).** If the deviation has fixed shape of amplitude τ^k over a fixed angular window (including all ends where the edge normals are exponentially pinned), the coefficients are period-type constants (ln 2, ln(1+√2), λ-logarithms, Catalan-type integrals) — no zeta values at any order.
**R3 (parity → absence of ζ(3)).** ζ(3) needs the asymmetric nesting: the tent's |s|-slope times the profile, or the profile's exponential square. A symmetric/rigid/even profile (log cosh, 1/cosh, flat Minkowski shift) has zero odd moments: ζ(3) never appears — even with algebraic overshoot (soft-max square) and even with a genuine layer.

**Test matrix (all confirmed):**

| family/end | layer profile | prediction | constants found |
|---|---|---|---|
| squircle p→1 | softplus, shrinking, overshoot | R1: zeta at every order | ζ(2)τ², ζ(3)τ³, ζ(4)τ⁴, ζ(5)τ⁵, ζ-products τ⁵–τ⁷, MZV at τ⁸ ✓ |
| squircle p→∞ | non-shrinking (exp-narrow layer) | R2 | period constants, no ζ ✓ |
| Lamé λ=2,3,1/2, p→1 | exp-narrow physical layer, analytic bulk | R2 | analytic constants, no ζ ✓ |
| Lamé p→∞ | non-shrinking F(θ) (closed form) | R2 | elementary d₁, no ζ ✓ |
| soft-max square | even bump (√2/p)F(β), pinned corners | R2+R3 | 16 ln2 − 8√2 ln(1+√2), single 1/p term, no ζ ✓ |
| LSE diamond | same, rotated | R2+R3 | exact single 1/p term, no ζ ✓ |
| hexagon p-norm | pinned, non-shrinking bump | R2 | period constants, no ζ ✓ |
| octagon p-norm | pinned, non-shrinking bump | R2 | period constants, no ζ ✓ |
| Minkowski sum | flat (no layer) | R3 | 2πr only, no ζ ✓ |
| sech family | rigid even profile | R3 | 4πw² only, no ζ ✓ |

---

## Files (work dir: new_results/AV5_layer_universality)

- `progress.log` — step log.
- `step_01.py … step_117.py` (harness-numbered). Key steps: 03–10 squircle pipeline & fits; 11–17 layer expansion (sympy g's, verified numerically); 18–27 f_k monomial assembly; 19–37 I(a,b)/K integrals (dps 90/135/210, tail-fix at s=80); 38–41, 47–50 PSLQ identification with pigeonhole bounds and two-precision residual verification; 51–52 high-precision F₇/F₈ + zeta-collapse boundary; 54–55 exact area; 56–57 square end; 58–79 anisotropic Lamé (support-function method, profile closed forms, e^{σ}/fan tests); 80–85 soft-max square & LSE diamond (exact h, exact c₁); 86–99 hexagon (Legendre-dual support function, NaN-poisoning fix, profile extraction); 100–117 octagon (45°-symmetric polar form, coefficients, pinning), family 4.
- Checkpoints: `ckpt_sq2_taus.json/ckpt_sq2_Ps.json` (squircle perimeter data), `ckpt_monomials.pkl`, `ckpt_integrals135.pkl`, `ckpt_integrals210.pkl`, `ckpt_Fk_125.pkl`, `ckpt_Fk78_200.pkl`, `pslq_I_120.json`, `ckpt_lame_p1.json`, `ckpt_lame_p1_coeffs.json`, `ckpt_lame_pinf_coeffs.json`, `ckpt_softmax.json`, `ckpt_hex3_*.json`, `ckpt_oct4_*.json`. (Superseded buggy data: `ckpt_hex_Ps.json`, `ckpt_hex2_*`, early `ckpt_oct*` — from the NaN-poisoned run.)

## Open items
- Closed forms for c₈, c₉ of the squircle diamond end (contain weight-8 alternating-MZV irreducibles; values given to ~58 digits, collapse boundary established at dps 200).
- Closed forms for the non-zeta period constants (squircle square-end d₁…, hexagon/octagon d₁…, anisotropic p→1 c₁…).
- Hexagon/octagon area expansions (cheap by polar integration; not completed).
- Anisotropic p→1 fan-boundary shoulder: detailed multiscale profile not fully resolved (expansion established as analytic).