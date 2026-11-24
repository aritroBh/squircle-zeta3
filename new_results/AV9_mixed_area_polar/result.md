# AV9 — The mixed area V(K_p, K_p′) of the squircle and its Hölder polar: exact duality in place of parity

**Object.** K_p = {(x,y): |x|^p+|y|^p ≤ 1} (the squircle), p′ = p/(p−1) its Hölder conjugate, K_p° = K_{p′} its polar. We study V(p) = V(K_p, K_{p′}) (2-d mixed volume/mixed area) at the square end (η = q−1 → 0, q = p′) and at the circle end (s = 1/p − 1/2 → 0), plus the broken-symmetry family V(K_p, K_r), r ≠ p′.

**Master formulas (PROVED from support-function theory; verified numerically).** With q = p′, u = r′,
- V(K_p,K_r) = (2/r) ∫₀¹ [v(1−v)]^{−1/u} (v^{q/u}+(1−v)^{q/u})^{1/q} dv (boundary-of-K_r form);
- V(K_p,K_r) = 2(q−1) ∫₀^{π/2} (c^u+s^u)^{1/u} (cs)^{q−2} (c^q+s^q)^{1/q−2} dθ,
using the exact radius-of-curvature lemma: for h(θ) = (cos^qθ+sin^qθ)^{1/q},
**h + h″ = (q−1)(cosθ sinθ)^{q−2} (cos^qθ+sin^qθ)^{1/q−2}** (derived by hand: A+B′ = (q−1)c²s²(c^{q−4}+s^{q−4}), and (c^{q−4}+s^{q−4})(c^q+s^q) − (c^{q−2}−s^{q−2})² = (cs)^{q−4}(c²+s²)²; verified to 51 digits, step_02b).
The given self-dual formula V(p) = (2/q)∫₀¹(v(1−v))^{1/q−1}(v^{q−1}+(1−v)^{q−1})^{1/q}dv was re-derived and matches the support-function mixed area to 27–146 digits at p = 1.5, 2, 2.5, 3, 4, 5 (step_01). **Duality: V(η) = V(1/η) exactly** (mixed area is symmetric in its arguments; η ↔ 1/η is p ↔ p′). Verified to 50–69 digits at η = 1/20, 1/10, 1/4, 1/2, 1, 3 with a well-conditioned θ-form (fix5). V(2) = π, V(1) = V(∞) = 4.

---

## Q1 — c₁…c₄ verified; exact derivation; this is NOT a Gamma quantity

Expanding the log-integrand in η (cumulant series of ln(e^{ηℓ}+e^{ηℓ′}), ℓ = ln v, ℓ′ = ln(1−v)) and integrating termwise with the **exact Beta derivatives**
∫₀¹ ℓ^a ℓ′^b dv = a!b![x^a y^b] B(x+1,y+1), where ln B = Σ_{k≥2}(−1)^k ζ(k)/k·(x^k+y^k−(x+y)^k) − ln(1+x+y) (γ cancels identically), gives c_n ∈ ℚ[ln 2, ζ(2), ζ(3), ζ(5), …], total weight ≤ n. Computed symbolically (av9core.py, step_03.py; 6 s to order 8, 20 s to order 12):

c₁ = −4 ln 2 ✓, c₂ = 2 ln²2 + 4 ln 2 ✓, c₃ = π²/6 − 4 ln2 − 4 ln²2 − 2ln³2/3 ✓, c₄ = ζ(3) − π²/3 − π²ln2/6 + ln⁴2/6 + 2ln³2 + 4ln2 + 6ln²2 ✓ — all four plan coefficients match **exactly** (π²/6 = ζ(2), π²ln2/6 = ln2·ζ(2)).

**Status: PROVED (exact derivation) + NUMERICAL corroboration:** numerical quadrature of the coefficient functions agrees to 87–102 digits (step_04/15/17); Cauchy–FFT of V(η) on |η| = 0.01–0.02 agrees to 24–44 digits.

**Alternating-MZV, not Gamma (explicit statement, as requested):** the η-expansion produces integrals ∫ ln^a v ln^b(1−v) dv (Beta derivatives at integer points — Riemann ζ values) with powers of ln 2 coming from the symmetric sum v^η+(1−v)^η (the constant of its cumulant generating function). Hence c_n lies in the subring ℚ[ln 2, ζ(2), ζ(3), ζ(5), ζ(7), ζ(9), ζ(11), …] of the level-2 (alternating) MZV algebra. Since ζ(3), ζ(5), ζ(7) genuinely occur (c₄, c₆, c₈), **V(η) is not expressible in Γ-values at rational arguments** (those generate only algebraic multiples of π-powers). Moreover **no Li_k(1/2), k ≥ 4, ever appears** through weight 12 — the only alternating generator is ln 2. The "incomplete-basis failure at order 5" noted by the probe is therefore explained: the missing basis elements were ordinary ζ-products (π⁴, π²ln2², ln2·ζ(3) in c₅), not alternating polylogs.

## Q2 — Orders 5–12 and the first-appearance map

Exact c₅…c₈ (full forms in step_03 output and the table below; c₉…c₁₂ in scratchpad pickles, step_14/16, verified step_15/17). First appearances on the self-dual line:

| constant | first order | coefficient | status |
|---|---|---|---|
| ln 2 | 1 | −4 | PROVED |
| ζ(2) (as π²/6) | 3 | 1 | PROVED |
| ζ(3) | **4** | 1 | PROVED |
| ζ(2)ζ(3) (as π²ζ(3)) | 6 | −1/2 | PROVED |
| ζ(5) | **6** | 1/2 | PROVED |
| ζ(3)² | 7 | −1/4 | PROVED |
| ζ(7) | **8** | 1/2 | PROVED |
| ζ(3)³ | 10 | 1/12 | PROVED |
| ζ(3)ζ(5) | 10 | 11/8 | PROVED |
| ζ(9) | **10** | −1/12 | PROVED |
| ζ(11) | **12** | 11/4 | PROVED |

**ζ(2k+1) first appears at order 2k+2 (k = 1…5 verified; conjecture: all k).** ζ(n) itself never appears in c_n (proved below). c₁…c₉ alternate in sign, then c₁₀ ≈ −1.889 breaks alternation and |c₁₂| ≈ 50.6 grows — the Taylor series converges only in a disk |η| < ρ with ρ estimated 0.7–1 from |c_n|^{1/n} (OPEN: locate the nearest singularity; it is not η = −1).

The table (exact; numerics to 31 digits verified as above):

c₅ = −8ln²2 − 3ζ(3) − 4ln2 − 4ln³2 − ln2 ζ(3) − 2ln⁴2/3 − ln⁵2/30 + π⁴/720 + π²ln²2/12 + π²ln2/2 + π²/2 = −3.6611143915541052553202258682395
c₆ = ζ(5)/2 + 6ζ(3) + 4ln2 ζ(3) + ln²2 ζ(3)/2 − π²ζ(3)/12 + ln⁶2/180 + ln⁵2/6 + 5ln⁴2/3 + 20ln³2/3 + 10ln²2 + 4ln2 − π²ln2 − 2π²/3 − 29π⁴/1440 − π²ln²2/3 − π⁴ln2/720 − π²ln³2/36 = 3.4247806814020087389447841948572
c₇ (contains −ζ(3)²/4, −13ζ(5)/4, −10ζ(3), …, full form in step_03) = −2.5554933609523348336049371142189
c₈ (contains ζ(7)/2, 9ζ(3)²/8, 43ζ(5)/4, 15ζ(3), …) = 3.9417879517107780986898404031641
c₉ (no ζ(9); ζ(7) coeff −23/8 − ln2/2) = −5.8783008607478681413956026111238
c₁₀ (ζ(9) first: −1/12; ζ(3)³/12; 11ζ(3)ζ(5)/8) = −1.88868798658935180146851277649
c₁₁ = 7.372586619120823262601324 (no ζ(11) ✓), c₁₂ = 50.626165468448716764202677797 (ζ(11) first: 11/4).

## Q3 — Why no ζ(3) at order η³: the top-weight theorem (PROVED); it is NOT the duality

**Theorem (proved, all n; verified n ≤ 12).** For the self-dual square-end expansion (and in fact for any family P(η)∫₀¹(v(1−v))^{α(η)}(v^{β(η)}+(1−v)^{β(η)})^{γ(η)}dv analytic at 0 with α(0) = β(0) = 0):
(i) every term of the coefficient functions h_k(η) has (ℓ-degree) + (ln2-degree) ≤ k, and the "extended top part" of e_n = [η^n]exp(Σ h_k η^k) equals
T_n = (−1)^n [(ln 2 + ℓ)^n + (ln 2 + ℓ′)^n] / (2 n!);
(ii) T_n contains only **pure** powers of ℓ, ℓ′, whose integrals ∫ℓ^j = (−1)^j j! are rational — so the entire weight-n zeta-part of c_n vanishes, in particular **[ζ(n)] c_n = 0 for all n**;
(iii) the weight-n part of c_n is exactly **4 (−ln 2)^n / n!** (i.e. the symbol of V(η) in the weight-graded MZV algebra is 4·2^{−η}).

*Proof.* h_k = Σ_{j≤k}(−1)^{k−j}g_j + (−1)^k(ℓ+ℓ′) + (−1)^k/k with g_j the cumulants of (e^{ηℓ}+e^{ηℓ′})/2; the weight-1 part of h₁ is −ln2 − (ℓ+ℓ′)/2 and the weight-k part of h_k (k ≥ 2) is g_k, so T_n = [η^n] exp(−ln2·η − (ℓ+ℓ′)η + Σ_{k≥2} g_k η^k) = [η^n] e^{−η ln2}(e^{−ηℓ}+e^{−ηℓ′})/2, using Σ_{k≥1} g_k η^k = ln((e^{ηℓ}+e^{ηℓ′})/2). All weight-n contributions come from T_n (ζ-weight ≤ ℓ-degree, so weight n forces ℓdeg+ln2deg = n). Verified: the extended top part of e_n equals T_n for n = 1…8 (verify_Tn.py). ∎

**Consequence:** ζ(3) (weight 3) cannot occur in c₃ — its absence is *forced by the integrand structure (single power-sum + Beta(1,1) measure), not by the p ↔ p′ duality*. Confirmed explicitly: in the broken one-parameter families W_c(η): u = 1+c/η (step_05) and r = q^c (step_06), d₁, d₂ are c-independent, and **[ζ(3)] d₃(c) = 0 identically**, [ζ(5)] d₅(c) = 0 identically, with d_n(1) = c_n for all n ≤ 6 (exact consistency check passed). So at the square end, breaking the duality does **not** lower ζ(3)'s first appearance below 4. The duality's genuine ζ(3)-suppressing action is at the circle end (Q4/Q5).

## Q4 — Circle end: V(s) = π + Σ a_{2k} s^{2k}; no ζ(3) at all through s¹⁰

With q = 2/(1−2s), exact termwise integration uses half-integer Beta derivatives (ψ^{(m)}(1/2) = (−1)^{m+1}m!(2^{m+1}−1)ζ(m+1), ψ(1/2) = −γ−2ln2; γ cancels identically; step_07b/10). Results:

a₁ = a₃ = a₅ = a₇ = a₉ = **0 exactly** (forced by V(η) = V(1/η), i.e. V even in s — PROVED);
a₂ = 3π(π² − 8)/4 = 4.4051515886861057008313760006490
a₄ = π(15π⁴ − 352π² + 1664)/192 = −5.7099163755880153855453367874535
a₆ = π(9699π⁶ − 96280π⁴ + 13440π² − 39936)/11520 = 10.547784674451002002666058478381
a₈ = π(−418761675π⁸ + 4482752384π⁶ − 3451302400π⁴ − 63651840π² + 138117120)/77414400 = −17.927818329122476877090216238499
a₁₀ = π(434288771763π¹⁰ − 4895069752920π⁸ + 6167864668672π⁶ − 1569660825600π⁴ − 5171281920π² + 11312824320)/6502809600 = 26.001119867216311827293139211912

**Every coefficient through s¹⁰ lies in π·ℚ[π²]: no ζ(3), ζ(5), ζ(7), ln 2, G anywhere** (raw expressions contain thousands of ψ^{(m)}(half-integer) terms with ζ(3), ζ(5), ln 2 — they all cancel). Verified by Cauchy–FFT of V(s) to 38–48 digits (step_08, and a₁₀ matches FFT exactly). This is far stronger than the probe's expectation ("ζ(3) at s⁶ or s⁸ if at all"): **the plan's Q4 answer is that ζ(3) never enters through s¹⁰; conjecture (NUMERICAL): V(s)/π ∈ ℚ[[s², π²]] to all orders.** Method note: exact Beta-derivative series replaced the 10-minute order-6 numerical attempt; a₀…a₁₀ took ~6 min total.

## Q5 — Broken symmetry at the circle end: ζ(3) enters at the first possible order

Bivariate expansion of V(K_{p₁}, K_{r₂}) with p₁ = 2/(1+2s₁), r₂ = 2/(1+2s₂) to total order 5 (step_09c; coefficients b_{jk} ∈ π·ℚ[ln2, ζ(3), ζ(5)]; verified by 2-D finite differences of the θ-integral to 11–17 digits, step_11; swap symmetry b_{jk} = b_{kj} exact; reproduces the self-dual line s₂ = −s₁ exactly, including a₂, a₄; reproduces the area line s₂ = s₁ = area(K_{p₁}) = 4Γ(1+1/p)²/Γ(1+2/p) exactly through order 4, ctrl.py):

- ln 2 enters at **order 1**: b_{1,0} = b_{0,1} = π(1 − 2 ln 2).
- ζ(3) enters at **total order 3**: b_{0,3} = b_{3,0} ⊇ −23πζ(3)/8 (this axis is the perimeter line: V(disk, K_r) = perimeter(K_r)/2), b_{1,2} = b_{2,1} ⊇ +15πζ(3)/8.
- ζ(5) enters at total order 5: b_{0,5} ⊇ −861πζ(5)/128.
- Along s₂ = λs₁, the order-3 ζ(3) coefficient is **−(λ+1)(23λ² − 38λ + 23)/8**; the quadratic has only the complex roots (19 ± 2i√42)/23, so **the self-dual line λ = −1 is the unique real direction with no order-3 ζ(3)** (the order-5 ζ(5) polynomial −3(λ+1)(287λ⁴−92λ³−262λ²−92λ+287)/128 likewise vanishes in real λ only at λ = −1). Odd orders vanish on λ = −1 by the swap symmetry; and the even orders are pure-π (Q4).

So breaking the Hölder duality lets ζ(3) in at order 3 (the lowest order compatible with the ℓ-degree bound at the circle end), while the exact duality suppresses it to all checked orders. This is the direct, quantitative confirmation of the paper's parity mechanism.

## Controls (required)

- **Parity-preserving control (must NOT produce ζ(3)): the self-dual circle-end expansion V(s).** It does not: odd orders vanish exactly (duality), and all even coefficients through s¹⁰ are pure π·ℚ[π²] — the precise analog of "Parseval gives only ζ(2n)". ✓ behaves as required.
- Secondary: square-end c₃ has no ζ(3) (top-weight theorem; not a parity effect).
- **Positive controls (should produce ζ(3)): the broken directions.** Area line s₂ = s₁ (a Γ-quantity, computed only as a control, not re-derived as a result): order-3 coefficient contains +8ζ(3)/3·(…) ✓; perimeter axis s₁ = 0: −23πζ(3)/8 at order 3 ✓; generic direction λ ≠ −1: ζ(3) at order 3 ✓. Also, the square-end broken families correctly do NOT produce ζ(3) before order 4 (structural top-weight obstruction), matching the theorem.

## Conjectures (deliverable 7c)

1. (PROVED) V(η) = V(1/η); V(s) even; odd circle-end coefficients vanish.
2. (PROVED) Top-weight theorem: weight-n part of c_n is 4(−ln2)^n/n!; ζ(n) never in c_n.
3. (NUMERICAL, orders ≤ 12) ζ(2k+1) first appears in c_{2k+2} with coefficient sequence 1, 1/2, 1/2, −1/12, 11/4 (k = 1…5) — **OPEN: identify this sequence / the generating function of the sub-top weights.**
4. (NUMERICAL, through s¹⁰) Circle end: V(s)/π ∈ ℚ[[s², π²]] — only π and rationals, all orders. **OPEN: prove it, or find the closed form** (coefficients of π^{2j} as rational series in s² are tabulated in step_18; no elliptic/gamma closed form identified).
5. (NUMERICAL, orders 3, 5) at the circle end the self-dual line is the unique real direction killing the leading odd-zeta term.

## Paragraph for the paper (deliverable 7d)

"Sections VIII–XV asked whether a non-linear mechanism can let ζ(3) past the parity barrier that excludes it from the Fourier/Parseval computations. The mixed area V(K_p, K_{p′}) of the ℓ^p-ball and its polar makes this precise. The exact Hölder-duality V(K_p,K_{p′}) = V(K_{p′},K_p) — an equality of the two arguments, not an imposed parity — forces V to be even in s = 1/p − 1/2. At the symmetric point p = 2 this symmetry is exactly the parity operator of Sections II–VII, and its effect is total: V(2e^s-like) = π + Σa_{2k}s^{2k} with every coefficient in π·ℚ[π²] through s¹⁰ (and conjecturally always) — no ζ(3), no ln 2, nothing but the even-zeta content that Parseval also produces. The moment the duality is broken, V(K_p,K_r) with r ≠ p′, the suppression fails at the first opportunity: ζ(3) appears at total order 3 in every direction of the (p,r)-plane except the self-dual line, which is the unique real line on which the order-3 ζ(3)-polynomial −(λ+1)(23λ²−38λ+23)/8 vanishes. (At the square end p → ∞ a different, symmetry-independent mechanism — a top-weight cancellation in the Beta-integral expansion — delays ζ(3) only to order 4, showing that at non-symmetric points ζ(3) is present regardless: parity is what excludes it, and exact duality is what enforces parity.)"

## Evidence table

| result | method | check | digits |
|---|---|---|---|
| given V-formula = mixed area | derivation | support-function mixed area, 5 p-values | 27–146 |
| R(θ) = (q−1)(cs)^{q−2}(c^q+s^q)^{1/q−2} | derivation | numeric h+h″ | 51 |
| c₁…c₈ exact | symbolic (sympy) | quad of coefficient functions | 100–102 |
| 〃 | 〃 | Cauchy–FFT of V(η) | 24–32 |
| c₉…c₁₂ exact | symbolic | quad / FFT | 87–90 / 31–44 |
| V(η) = V(1/η) | mixed-area symmetry | 6 η-values | 50–69 |
| top-weight T_n identity | proof in text | explicit n ≤ 8 | exact |
| a₂…a₈ exact | symbolic | Cauchy–FFT of V(s) | 38–48 |
| a₁₀ exact | symbolic | FFT | exact match |
| odd circle coeffs = 0 | duality | a₁,a₃,a₅,a₇,a₉ symbolic + FFT | exact/1e-34–1e-50 |
| b_{jk} (order ≤ 5) exact | symbolic | 2-D finite differences (11 pts) | 11–17 |
| b_{jk} swap, self-dual & area lines | symbolic | exact identities | exact |
| broken families d_n(c) | symbolic | d_n(1) = c_n, n ≤ 6 | exact |

**PSLQ record (deliverable 8; mpmath pslq, basis = {ln2^a ζ(2)^b ζ(3)^c ζ(5)^d ζ(7)^e, weight ≤ n}):** c₅: exact match at dps 100 & 150, max coeff 480, resid 1.8e-99 / 0; c₆: match at dps 150 & 200, max 3600, resid 1.6e-148 / 5.0e-198; c₇: match at dps 300, max 30240, resid 7.7e-298; c₈: not completed within time cap (exact primitive max 201600; dps 250/300 runs returned overfit near-relations, resid ~1e-188, that fail the two-precision rule). **Cautionary finding (relevant to the plan's "incomplete-basis failure"):** at insufficient precision PSLQ returns spurious relations (c₆@dps100, c₇@dps100/200, c₈@dps250); all exact forms here are instead derived symbolically and independently verified to the digits above. Incomplete-basis controls (c₅ w/o ζ(3), c₆ w/o ζ(5), c₈ w/o ζ(7)) correctly return no sub-1e4 relation.

## Files (cwd = new_results/AV9_mixed_area_polar)

av9core.py (exact series/Beta machinery); step_01.py (formula vs mixed area); step_02b.py (R(θ), bivariate formula); step_03.py (c₁…c₈); step_04.py (quad+FFT verify); step_05.py, step_06.py (broken square-end families); step_07b.py (circle end a₀…a₈); step_08.py (circle FFT verify); step_09c.py (bivariate circle end + 3 line checks); step_10.py (a₁₀); step_11.py (FD verify); step_12c.py, pslq2.py, pslq3.py (PSLQ); step_13b.py/fix3.py (top-weight); fix5.py (duality, well-conditioned); step_14.py, step_16.py (c₉…c₁₂, first appearances); step_15.py, step_17.py (verification); step_18.py (summary tables); ctrl.py (area-line control); verify_Tn.py (T_n identity); scratchpad/*.pkl, pslq_c*.txt (checkpoints); progress.log (one line per step).

**Bottom line for the parity model.** The squircle mixed area gives the model a exact mathematical home: the Hölder duality *is* the parity. Where it holds (circle end), the expansion is purely even-zeta (π·ℚ[π²]) to all orders computed — ζ(3) is absent not by accident but by symmetry; where it is broken, ζ(3) enters at the first order the weight bound allows (order 3), and the self-dual line is the unique exempt direction. At the non-symmetric square end, ζ(3) is present from order 4 on (ζ(5) from 6, ζ(7) from 8, ζ(9) from 10, ζ(11) from 12) on every path, dual or not — its order-3 absence there is a separate, proved, symmetry-independent top-weight cancellation.