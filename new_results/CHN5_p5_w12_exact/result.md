# CHN5: an exact weight-12 cut of the P5 curve Γ (Kimi K3)

## Question
Does Γ⁺ (the open-orthant part of Γ = V10⁺ ∩ {P7 = 0}) contain a point with **A = CLASS = 0**?
(A = coef ζ2²ζ(3,5) of c₁₂; B = coef ζ2ζ(3,7); CLASS = coef ζ(3,9) + (14/9)·coef ζ(5,7); Z = coef ζ(2,1,1,8).)

## VERDICT: **NO** — A < 0 on *all* of Γ⁺, so no interior point is irreducible-free through weight 12.

- **EXACT (PROVED):** A, B, CLASS, Z each factor as a₂a₃a₄(Σa)·(degree-8 quotient); on V10 the W12 irreducible content is spanned by **A and CLASS alone** (B = (75/14)A, Z = −(25/3)A). Hence on Γ an interior point is irreducible-free ⟺ A = CLASS = 0.
- **NUMERICAL (very strong):** A = a₂a₃a₄(Σa)·Q_A < 0 and CLASS = a₂a₃a₄(Σa)·Q_C > 0 at **every** interior Γ point traced (a₁/a₃ from 10⁻¹¹ to 10⁹⁺, all branches found). The first factor is positive in the open orthant, so **A < 0 and CLASS > 0 on all of Γ⁺** — the two never vanish, together or separately, in the interior.

So, as expected, **P5 has no irreducible-free exponent vector in the open positive orthant through weight 12.**

---

## 1. Exact W12 forms — **EXACT (PROVED)**

Each of the four W12 irreducible coordinates of c₁₂ is a₂a₃a₄(Σa) times an irreducible reversal-invariant degree-8 polynomial:

| coord | formula | terms | coeffs +/− | P8(1,1,1,1,1) |
|---|---|---|---|---|
| A (ζ2²ζ(3,5)) | (1/100)·a₂a₃a₄(Σa)·**P8_A** | 367 | +243/−124 | 65888 |
| B (ζ2ζ(3,7)) | (1/224)·a₂a₃a₄(Σa)·**P8_B** | 359 | +187/−172 | 105168 |
| CLASS (ζ(3,9)+(14/9)ζ(5,7)) | (1/432)·a₂a₃a₄(Σa)·**P8_CLASS** | 450 | +223/−227 | 98600 |
| Z (ζ(2,1,1,8)) | (1/3)·a₂a₃a₄(Σa)·**P8_Z** | 226 | +154/−72 | 3712 |

Each P8 is integer-primitive, **irreducible over ℚ** (`factor_list`: single factor) and **reversal-symmetric** (monomial check, `subs(simultaneous=True)`). Files `P8_{A,B,CLASS,Z}.txt`.

**Proof of the factorisations.** The coordinates vanish exactly on each hyperplane a₂=0, a₃=0, a₄=0, Σa=0 (4 generic pts each), so a₂a₃a₄(Σa) divides them. The degree-8 quotients were interpolated on a rank-255 reversal-invariant determining set (255 orbits), then **upgraded to a full proof**: D = coord − a₂a₃a₄(Σa)·Q8 vanishes at **all 924 points** of a rank-924 determining set for reversal-invariant degree-12 polynomials (rank verified = 924 mod two primes, 2⁶¹−1 and 2⁶³−25, via `nmod_mat`). The m=0 (γ-prefactor) contribution to A, B, CLASS, Z is **exactly 0 at all 924 points**. All four identities independently re-verified against the engine at 3 fresh rational points (exact match). ⇒ identities **PROVED**.

## 2. The planner's proportionalities on V10 — **EXACT (PROVED)**

By exact ideal membership (linear algebra over ℚ, `linalg_tools.ideal_member`, verified):
- Q_B − (75/14)·Q_A = U₁·C4 + V₁·Q7′ (U₁ deg 4, 50 terms; V₁ deg 2, 6 terms)
- Q_Z + (25/3)·Q_A = U₂·C4 (U₂ deg 4, 50 terms; V₂ = 0)

Since a₂a₃a₄(Σa) > 0 in the orthant, on V10 = {C4 = Q7′ = 0}: **B = (75/14)·A and Z = −(25/3)·A exactly.** (Confirmed numerically to ~150 digits at refined Γ points.) These are *not* polynomial identities (they fail at generic rational points — planner_check.out in CHN4). ⇒ On V10 the W12 irreducible content = span(A, CLASS); an interior Γ point is irreducible-free ⟺ A = CLASS = 0.

(LP note: *no* sign-definite representative of Q_A or Q_C mod (C4,Q7′) exists on the *whole* open orthant — LP infeasible, formulation validated — so the sign must be argued on V10 itself, as done below.)

## 3. A on Γ⁺ — **NUMERICAL (very strong)**

Γ = {C4 = Q7′ = P7 = 0} (a₃ = 1 chart), traced by V10-sheet Q11-sign-change brackets + Newton on the exact polynomials (mpf, 60–150 digits; the trace uses only the *proved* exact forms, so evaluation is cheap and reliable).

**Result: A < 0 and CLASS > 0 at every interior Γ point, with no exceptions.**

- **119 interior points**, a₁/a₃ ∈ [10⁻⁷, 10⁶] (main sweep), branch counts 1–5 per a₁; **all** have sign(Q_A) = −, sign(Q_C) = +. Normalised min |Q_A|/S⁸ = 8.2·10⁻¹³ (at the a₁→0 end); in the interior it is ~10⁻⁷–1.9·10⁻⁶, **never 0**.
- **Small-a₁ end:** branch → stratum point (0, 8/3, 1, 8/3, 1); Q_A/S⁸ ≈ **−8.16201·10⁻⁶·a₁** (exactly linear over a₁ = 10⁻⁷…10⁻¹¹, negative slope), Q_C/S⁸ ≈ **+1.606·10⁻⁴·a₁**. Clean approach, no sign crossing; Q_A = 0 only *at* the boundary stratum.
- **Large-a₁ end:** Γ persists to a₁ ≥ 10⁹ (approaching the a₃→0 face); Q_A/S⁸ → negative plateaus ≈ −1.89·10⁻⁶ and −2.33·10⁻⁵. (The a₃=0 face kills A via the a₂a₃a₄Σa factor, not via Q_A.)
- **Stratum artifacts excluded:** at a₁ = 1 exactly, the boundary stratum {a₅=0, a₂=a₄, a₁=a₃} is a solution *curve* (Q_A = Q_C = 0 on it); 4 Newton landings on it (a₂=a₄ to 60 digits, a₅ ~ 10⁻⁶¹) were removed. The genuine nearby interior branch has Q_A < 0.
- **Cross-checks:** 11 independent dense solves (a₁ = 10⁻⁵…3·10⁶) all agree (Q_A<0, Q_C>0).

**Resultant corroboration (EXACT structure):** Res_{a₅}(C4,P7) = a₁⁴(a₁+1)⁴·F₃·F₃·F₈·F₆; Res_{a₅}(C4,Q_A) = a₁(a₁+1)·G30 (G30 deg 30, 4243 terms). Q_A = 0 at a Γ point forces G30 = 0 at its (a₁,a₂,a₄)-projection; **G30 ≠ 0 at all 119 traced points** (min |G30| = 6·10⁴). Traced points lie on G22 (residuals ~10⁻¹²⁸) and on F₆ (small-a₁ branch) / F₈ (larger branches).

**Exact-form validation:** the polynomial A evaluated at refined Γ points reproduces the independent mpf engine exactly (e.g. A = −46.7202620046 at a₁ = 0.01 vs engine −46.720262; CLASS = 919.264414442; B/A = 75/14 and Z/A = −25/3 to 150 digits).

**Remaining gap (why not PROVED):** a full exact elimination of {C4, Q7′, P7, P8_A} is infeasible here — the a₁-projection has degree ~10³ (Bézout 4·6·7·8 = 1344, plus two boundary-stratum curve components), and iterated resultants inflate to univariate degree ~10⁵–10⁶ (Sylvester 40×40+ with polynomial entries). The step-3 sign statement is therefore labelled NUMERICAL, but it is multi-method, covers all branches found over 20 decades in a₁/a₃, has clean asymptotics at both ends, and is corroborated by the exact resultant condition G30 ≠ 0.

## Evidence table

| check | method | result | residual/digits |
|---|---|---|---|
| W12 basis indices (A=7,B=8,CLASS=9+14/9·10,Z=11) | mzv_exact BASIS | confirmed | exact |
| reversal-invariance of all 12 coords | 2 generic Fraction pts | invariant | exact |
| m0 contribution to A,B,CLASS,Z | 924 det pts | **all exactly 0** | exact |
| hyperplane vanishing (a₂,a₃,a₄,Σa = 0) | 4 pts each, full engine | all 4 coords = 0 | exact |
| det12 rank (reversal-inv deg-12) | nmod_mat mod 2 primes | **924** | exact |
| A,B,CLASS,Z = a₂a₃a₄(Σa)·Q8 | D=0 at all 924 pts | **PROVED, 0/924** | exact |
| quotient irreducibility + reversal-symmetry | factor_list / monomials | irreducible, symmetric | exact |
| fresh-pt re-verification (3 random pts) | exact engine | all 4 match | exact |
| B−(75/14)A ∈ (C4,Q7′) | ideal membership | U₁C4+V₁Q7′, verified | exact |
| Z+(25/3)A ∈ (C4,Q7′) | ideal membership | U₂C4 (V₂=0), verified | exact |
| Q_A, Q_C signs on Γ⁺ (119 pts, a₁/a₃ 10⁻⁷–10⁹) | mpf 60–150d trace | **all Q_A<0, Q_C>0** | residuals ≤10⁻¹³³ |
| small-a₁ asymptotic | exact polys, a₁=10⁻⁷–10⁻¹¹ | Q_A/S⁸~−8.16201e-6·a₁ | 60 digits |
| exact A vs engine at Γ pts | 150-digit eval | matches (e.g. −46.720262) | 150 digits |
| G30 (Q_A-resultant) ≠ 0 on Γ⁺ | 119 traced pts | min |G30| = 6e4 | exact |
| stratum artifacts (a₁=1) | a₂=a₄, a₅~10⁻⁶¹ | excluded (Q_A=Q_C=0 on stratum) | — |

## Files and timings (work dir `CHN5_p5_w12_exact/`)

- **Exact forms:** `P8_{A,B,CLASS,Z}.txt`, `Q8_quotients.pkl`, `exact_forms.log`, `factor.log`; proof data `det12_points.pkl`, `eval12_*.pkl` (924 pts), `interp12.log` (0/924 for all four), `run_eval12.py`, `run_interp12.py`.
- **Proportionalities:** `ideal.log`, `ideal_B-75_14A.pkl`, `ideal_Z+25_3A.pkl`, `run_ideal.py`.
- **Γ analysis:** `gamma_trace2.pkl/.log` (mpf, a₁ 10⁻⁷–10⁶), `gamma_trace.pkl/.log` (first pass), `refined_gamma_pts.pkl`, `run_gamma_trace2.py`; resultants `resultant_{P7,QA}.pkl`, `resultants.log`, `run_resultants.py`.
- **Timings:** det-set build ~3 s; W12 exact eval ~7.2 s/pt (924 pts ≈ 5 chunks × ~1200 s ≈ 1.7 hr); interpolation+924-pt proof ~2 min; ideal membership <1 s each; quotient factorisation <1 s; Γ traces ≈ 55 min (numpy pass) + 33 min (mpf pass) + targeted probes; Res_{a₅}(C4,P7) 14 s, Res_{a₅}(C4,Q_A) 143 s.
- **Auto-saved step cells:** `step_01.py … step_115.py` (setup, det-set, campaign, interpolation, tracing, factor identification, validations).

## Labels
- **EXACT (computer algebra):** the four W12 factorisations (PROVED, 924-pt rank-924); quotient irreducibility + reversal-symmetry; the two proportionalities on V10 (ideal membership); V10 = {C4,Q7′} and f11 = −(1/10)a₂a₃a₄(Σa)P7 (from CHN4); the resultant factorisations Res_{a₅}(C4,P7), Res_{a₅}(C4,Q_A).
- **NUMERICAL (mpmath, 60–150 digits, multi-method with cross-checks):** Q_A < 0 and Q_C > 0 at every interior point of Γ⁺ (a₁/a₃ ∈ 10⁻¹¹–10⁹⁺), hence **no point of Γ⁺ has A = CLASS = 0**; Γ⁺ nonempty; stratum-artifact exclusion; both end-asymptotics.
- **CONJECTURAL (MZV independence):** needed only to read "coefficient ≠ 0" as "the irreducible is genuinely present in the volume"; all vanishing/nonvanishing coefficient statements themselves are exact.

**Bottom line.** Through weight 12, the only irreducible-free points of P5 in the closed orthant are on the boundary (faces a₂,a₃,a₄ = 0 and the strata {a₁=0,a₂=a₄,a₅=a₃}, {a₅=0,a₂=a₄,a₁=a₃}). On the W11-obstruction curve Γ = V10⁺∩{P7=0}, the weight-12 coordinate A is strictly negative everywhere in the open orthant (and CLASS strictly positive), so **Γ⁺ contains no point with A = CLASS = 0 — P5 has no irreducible-free exponent vector in the open orthant through weight 12** (answer: NO; polynomial content EXACT, on-Γ sign statement NUMERICAL but overwhelming).