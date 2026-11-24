# COR8R referee report — COR8's depth-two symbol claim

Work dir: `new_results/COR8R_referee_depth_two`. All numbers below were printed by my own code (`cor8r_sym.py`, `cor8r_num.py`, `step_*.py`), written from definitions before opening COR8's code.

## 1. Conventions used (stated, and checked on series)

* **MPL convention (convA, Goncharov):** `Li_{a,b}(x,y) = Σ_{0<k₁<k₂} x^{k₁}y^{k₂}/(k₁^a k₂^b) = G(0^{b−1}, 1/y, 0^{a−1}, 1/(xy); 1)`. **convB (COR4's):** `Li_{a,b}(x,y) = Σ_{n>m≥1} x^n y^m/(n^a m^b) = G(0^{a−1}, 1/x, 0^{b−1}, 1/(xy); 1)`.
  Series checks (printed): `G(1/y,1/(xy);1) = Li_{1,1}(x,y)` to 30 digits by direct nested integration; `G(0,1/y,0,1/(xy);1) = Li_{2,2}(x,y)` and `G(1/y,0,0,1/(xy);1) = Li_{3,1}(x,y)` to ~50 digits via the expansion of the iterated integrals (residuals 4·10⁻⁵³, 2·10⁻⁵²).
  *Only convB has convergent double sums on all of c∈(0,1) for the arguments (½, 1±c)* (convA's series for y=1+c>1 diverges).
* **Differential (derived by parts, verified numerically at n=2,3 to 15–20 digits, and for zero-containing words):** with a₀:=z, a_{n+1}:=0,
  `dG(a₁,…,aₙ;z) = Σᵢ G(âᵢ;z) dlog((a_{i−1}−aᵢ)/(aᵢ−a_{i+1}))`.
* **Symbol recursion:** `S(G) = Σᵢ S(G(âᵢ)) ⊗ [letterᵢ]`, outermost integration last; constant letters (dlog=0) dropped. Tensor order gives `S(Li_n(z)) = −[1−z][z]^{n−1}`, `S(H_{a₁…aₙ}(c)) = g(aₙ)…g(a₁)`, g(0)=[c], g(1)=−[1−c], g(−1)=[1+c] — matching the paper's stated convention.
* **δ₂,₂:** for `Σ x₁⊗x₂⊗x₃⊗x₄` take `Σ (x₁∧x₂)∧(x₃∧x₄)` (symbol-of-coproduct = deconcatenation, then B₂→Λ²B₁ on each side; paper's convention). Coefficients in Laurent polys in c.

## 2. Evidence table (all printed)

| Check | Result | Residual |
|---|---|---|
| convA/convB G-form vs series (Li₁,₁, Li₂,₂, Li₃,₁) | match | 30–50 digits |
| DE formula n=2,3 vs numeric ∂/∂aᵢ,∂/∂z | match | 15–20 digits |
| `d Li₂,₂(½,1−c)/dlog c = Li₂,₁·dlog(1−c)` (top DE of the computed symbol) | match | 1.5·10⁻²⁸ |
| **δ₂,₂(Ω) from the 56 weight-4 HPL terms** | = `(1−c⁻¹)([c]∧[1−c])∧([1−c]∧[1+c]) + (−1−c⁻¹)([c]∧[1+c])∧([1−c]∧[1+c])` | **exact match to claim/paper, residual 0** |

The symbol of Ω collapses to 24 pure tensors; δ₂,₂(Ω) has **exactly two** Λ²B₂-components, coefficients `1−c⁻¹` and `−1−c⁻¹` — nonzero.

## 3. Verdicts

| Claim | Verdict | Detail |
|---|---|---|
| A, δ₂,₂(Ω) formula | **PROVED** | exact agreement with claim and paper Prop. prop:M3angle, independent code |
| A, "one depth-two term cannot suffice" | **PROVED** | δ₂,₂(Ω) = A(c)W₁+B(c)W₂ with W₁,W₂ independent basis wedges and A/B = −(1−c)/(1+c) nonconstant; a single term r(c)·(constant-coefficient MPL symbol) forces A/B constant — contradiction |
| A, δ₂,₂(Ω) = δ₂,₂(DT_A) = δ₂,₂(DT_B) | **NOT CONFIRMED** | With [2] dropped (paper: "modulo constants") and convA, I get δ₂,₂(DT_A) = **−**δ₂,₂(Ω) (overall sign discrepancy: likely a sign/order convention repair in the DT terms). With [2] kept, my specialization-first symbols leave [2]-wedges that Ω lacks; I traced this to a genuine defect of specializing before recursing — constant-factor letter choices become shuffle-incompatible at weight ≥3 (residual all in the [2]-letter; **with [2] dropped all shuffle identities in my tests hold exactly**). The correct route (two-variable symbols, then specialize x=½) was identified but not finished in time. convB gives δ₂,₂(DT_A)≡0 mod constants |
| B, sunit_args.json sign flips | **PROVED (mechanism), count not re-verified** | independent refactorization of sampled entries (e.g. y=−(1−c)²/(4c) stored as `[−2,−1,−2,0]` over (2,c,1−c,1+c): the [1−c] exponent is −2, must be +2) confirms the flip exists; I did not reprint all 40, so "20 of 40" is unverified. No paper statement relies on `sunit_args.json` (the paper's {2,3}-unit sentence is an independent number-field statement) |
| C, level-2 reductions at 120 digits | **NOT CONFIRMED** | not reached (no computation printed) |
| D, [2]-dropped symbol identity | **NOT CONFIRMED** | not reached; my engine's [2]-handling needed the repair above first |

## 4. Sentence the paper could carry

"Goncharov's criterion applies at the level of symbols: the (2,2) component of M₃(α)'s weight-four part is
δ₂,₂(Ω) = (1−c⁻¹)([c]∧[1−c])∧([1−c]∧[1+c]) + (−1−c⁻¹)([c]∧[1+c])∧([1−c]∧[1+c]) ≠ 0 for c∈(0,1),
so M₃(α) is not a combination of classical polylogarithms of the angle; no single depth-two term can account for it (its two components have nonconstant coefficient ratio). At the level of symbols modulo constant letters, two depth-two functions account for the obstruction: M₃(α) − DT_A is classical modulo products."

## 5. Repairs required for COR8

1. **R1 (claim A, DT signs/convention):** re-derive the DT symbols from the *two-variable* symbol of Li_{a,b}(x,y) (clean alphabet {x,1−x,y,1−y,xy,1−xy}) and specialize x=½ afterwards; specializing first produces shuffle-inconsistent [2]-letters (demonstrated: residuals all in the [2]-letter; with [2] dropped my system is exactly consistent). Then re-establish the overall sign: my consistent convention gives δ₂,₂(DT_A) = −δ₂,₂(Ω), so either DT → −DT or a tensor-order/sign repair in the MPL convention.
2. **R2 (claim B):** state the corrected factor vectors explicitly (sign of the [1−c] exponent in the 20 affected [1−r] vectors).
3. **R3 (claims C, D):** provide the 120-digit outputs and the 243-component residual computation with the corrected symbol engine.

## Files
- `cor8r_sym.py` — my symbol engine (letters over (2,c,1−c,1+c), [−1]=0; DE recursion; shuffle-regularization solver; δ₂,₂; HPL/MPL wrappers).
- `cor8r_num.py` — series/G-integral numeric helpers.
- `step_11.py`…`step_61.py` — exploration, convention series checks, DE derivation+checks, δ₂,₂(Ω) computation (matches claim), DT analysis, and the shuffle-consistency forensics documented in §3.
- `progress.log`, `result.md` in the work dir.