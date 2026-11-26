# Lean verification checklist

Started 2026-11-25. Paper: draft 58. Lean 4 with Mathlib.

Status words: TODO, STATED (Lean statement written and audited against the paper), PROVED (builds, no `sorry`,
only the standard axioms), BLOCKED (says why), OUT (cannot be done in Lean today; stays on the existing checks).

## What Lean can and cannot check here

Lean checks a proof line by line down to the axioms. It can only do that for mathematics that exists in its
library, Mathlib. Mathlib has `π`, `log`, `Γ`, `ζ(2)=π²/6`, `ζ(4)=π⁴/90`, integrals and power series. It does
not have polylogarithms, multiple zeta values, harmonic polylogarithms, symbols, or interval arithmetic at
300 digits. Most of the paper's main theorems are closed forms in exactly those objects. So:

- Wave 1 checks every exact, finite claim that Mathlib can express now.
- Wave 2 builds a small dilogarithm library and checks `M₁` and `M₁(α)` with it.
- The other main theorems (M₂ to M₄, pure weight, radius zero, decay law, hexagon theorems) stay verified the way
  they are now: written proofs plus independent computer re-checks. Lean does not replace that.

A Lean proof only says the Lean statement is true. Whether the Lean statement is the paper's statement is a
human-style reading job, so every target has a separate audit step.

## Workflow (per target)

1. **State** (Opus): write the Lean statement, quoting the paper line it formalises.
2. **Audit the statement** (a second Opus, sees only the paper and the Lean statement): same claim? same
   hypotheses? same constants? Reject if weaker or different.
3. **Prove** (Sonnet first; Opus if Sonnet fails twice): fill the proof.
4. **Gate** (mechanical, Haiku): `lake build` passes; no `sorry`, `admit`, `native_decide` or new `axiom`;
   `#print axioms` shows only `propext`, `Classical.choice`, `Quot.sound`.
5. **Record** (Haiku): fill the row below. A target that fails stays TODO or BLOCKED with the reason. Nothing is
   marked PROVED by hand.

If Lean finds that a statement in the paper is false as written, that goes to Aritro before anything else.

Rules: nothing in `paper/` is edited by this work. Lean files live in `lean/Zeta3/`. One `lake build` at a time.

## Wave 1: exact claims Mathlib can express now

Files are in `lean/Zeta3/Wave1/`. Status comes from `tools/gate.py check`, run by the planner on 2026-11-25
(results in `lean/audit/*.gate.json`). "PROVED n" means all n statements in the file are proved, the statement
block is the audited one, and the only axioms are `propext`, `Classical.choice`, `Quot.sound`.

| # | File | Paper result | What Lean proved | What it does not cover | Status |
|---|---|---|---|---|---|
| T01 | `T01_Rhombus` | Thm `thm:rhombus` | The curve `|x/A|+|y/B|=1` is the four sides; each side has length `r`; the perimeter is `4r`; the perimeter integral at `p=1` is `4r` | The first three derivatives of `P` at `p=1` | PROVED 5 |
| T02 | `T02_SeriesDiamond` | Prop `prop:an`, Thm `thm:pw` | The τ-series of Thm 8 through `τ⁵`, rewritten in `d`, has exactly the coefficients `a₀…a₅` of Prop 1 (with Mathlib's `ζ(2)`, `ζ(4)`); a `C⁵` function with the Thm 8 expansion has those Taylor coefficients | That the perimeter has either expansion; `a₆`, `a₇`; the weight claims | PROVED 3 |
| T03 | `T03_Omega` | Prop `prop:an` | `V₃=7`, `V₅=1` as true iterated derivatives; `ω₃=1`, `ω₅=-3/4` | `ω_n≠0` for all odd `n`; that `ω_n` is the `ζ(n)` coefficient | PROVED 4 |
| T04 | `T04_SmallFacts` | Thm `thm:pw`, `thm:oct`, `thm:zero` | `Y₂=1`, `Y₃=1/(2-p)`; `-(N²-1)((N-1)³-2)/(3N²)` is `1/4` at `N=2`, `-16/9` at `N=3`, never `0`; `ρ≥π/2` with equality iff `A=B` | The expansions themselves | PROVED 7 |
| T05 | `T05_Cm` | Thm `thm:Cm` | `C₂=6√2`, `C₃=-9√3/8`, `C_m<0` for all `m≥3`; the corners lie on `∂K₁`; `2m·Z_e=C_m` for every edge, from the geometry | That `C_m ζ(3)` is the `ζ(3)` part of the third derivative | PROVED 5 |
| T06 | `T06_PolyR1` | Thm `thm:poly` | `R₁=0` at `m=2`; `R₁=log 2` at `m=3`; `P(1)=2√3` and `P(1)R₁=2√3·log 2` for the hexagon | The identity `P'(1)=P(1)R₁`; "only for `m=2`" | PROVED 4 |
| T07 | `T07_PhiSeries` | Thm `thm:flat` | `Σ c^{2n}/(4n²-1) = 1/2-(1-c²)/(2c)·artanh c` for `0<|c|<1`, and for `c=cos(α/2)` | The `M_k` expansions | PROVED 2 |
| T08 | `T08_PolygonSums` | Prop `prop:edgelaw`, `prop:order5`, Thm `thm:Cm` | For the regular `2m`-gon, all `m≥2`: the edges lie on `∂K₁`; the edge data `L=2s²/c`, `h=s`, `ξ₀=ξ₁=2s²` follow from the geometry; `ΣZ_e=C_m`, `ΣZ_e^{(5)}=C_m^{(5)}`, `ΣZ_e^{(2,3)}=C_m^{(2,3)}` | That `Z_e` etc. are the zeta parts of the derivatives | PROVED 9 |
| T09 | `T09_OctSeries` | Thm `thm:oct` | The stated `d`-series through `d⁴` follows from the stated τ-series of `log((2-p)S/(4√3))` | That `S` has either expansion; `τ⁵`, `τ⁶` | PROVED 1 |
| T10 | `T10_PiP` | Thm `thm:pip` | The main-part integral is `4·2^{-1/p}`; it gives `N₀`, `N₁` and the `log` parts of `N₂`, `N₃`, with an explicit `O(p⁻⁵)` bound | `E(p)`, so the `π²` and `ζ(3)` parts | PROVED 5 |
| T11 | `T11_Decimals` | Thm `thm:xnyn` | All 10 printed decimals of `β₁` and of `β₀` | `d₃` (needs `ζ(3)` numerically); the paper prints no decimals for `b₃`, `c₃` | PROVED 2 |
| T12 | `T12_EdgeLaw` | Prop `prop:edgelaw` | From the geometry: `ΣZ_e=6𝓔` for every parallelogram and `-6𝓔/8` for every affine-regular hexagon, matching the two displays | The expansions of `P(1+d)` themselves | PROVED 4 |
| T13 | `T13_Lambda` | Thm `thm:P2` | `Λ_m=(P_{m-1}-P_{m+1})/(2m+1)` is a polynomial in `cos θ`; `Λ₁=½sin²θ`; Legendre values at 0; `Λ_m(π/2)=0` for even `m` and the explicit odd values; the Gamma closed form of the ε-term at `θ=π/2` (Bailey's summation, proved from scratch by a contiguous relation and Euler's limit formula) | The Mellin identity and the series for `L_θ(n)` | PROVED 9 |
| T14 | `T14_Rulers` | Prop `prop:ruler` | `L₁(p)=8` and `L_∞(p)=8·2^{-1/p}` for all `p≥1`, as true lengths (total variation) of the closed curve, and as the paper's integrals; `L_q→L_∞` as `q→∞` | Convergence of the `1/p` expansion of `L_q` | PROVED 6 |

Totals: 14 files, 66 statements, all 66 proved, 4649 lines of Lean. `lake build` of the whole library passes. No statement of the paper was found false.

## Wave 2: a dilogarithm library, and what follows from it

Files are in `lean/Zeta3/Wave2/`. Status from `tools/gate.py check`, run by the planner on 2026-11-25
(`lean/audit/*.gate.json`). The library defines the real dilogarithm as `Li₂(x) = -∫₀ˣ log(1-t)/t dt`.

| # | File | Paper result | What Lean proved | What it does not cover | Status |
|---|---|---|---|---|---|
| D00 | `D00_Dilog` | library (used in the proof of `lem:dilog`) | The integrand is integrable; `Li₂` is continuous on `(-∞,1]`; its derivative; the power series for `|x|≤1`; `Li₂(1)=π²/6`, `Li₂(-1)=-π²/12`, `Li₂(½)=π²/12-½log²2`; reflection, Landen, both inversions, duplication, Legendre's `χ₂` identity, Abel's five-term relation | The complex dilogarithm; `Li_n` for `n≥3` | PROVED 16 |
| D01 | `D01_TandM1Algebra` | Thm `thm:M2` eq. (T); end of the proof of `thm:M1` | `(π²+6λ²)/24=(ζ(2)-Li₂(½))/2` with Mathlib's `ζ(2)`; the closing algebra of the `M₁` proof gives the closed form of `thm:M1` | The six `I(·,·)` values are taken as printed | PROVED 5 |
| D02 | `D02_DilogSqrt2` | Lemma `lem:dilog` | All six identities at `b=√2-1`; every argument is `≤1` | | PROVED 7 |
| D03 | `D03_M1RightAngle` | Thm `thm:M1`, eq. `eq:M1t`, Prop `prop:M1short` at `α=π/2` | The integral `eq:M1t` is integrable and equals the closed form of `thm:M1`; the four-dilogarithm form at `α=π/2` is the same number; `8M₁` agrees with the printed `2.02095276965635835695` to `10⁻²⁰` | That `M₁=-F'(0)/2` is the integral `eq:M1t` (that is `lem:holo`, wave 3) | PROVED 4 |
| D04 | `D04_MalphaForms` | Thm `thm:Malpha` with Table `tab:corner`; Prop `prop:M1short` | For every `0<α<π`: `M₁(α)`, defined as the `u`-derivative at `0` of the paper's corner integral `M(u;α)`, equals the 19+8-dilogarithm table form and the four-dilogarithm form; every `Li₂` argument is `<1`; every logarithm argument is positive | The full expansion of `M(u;α)` in `u`; `M_k(α)` for `k≥2` | PROVED 5 |
| U01 | `U01_CubeIntegral` | Prop `prop:cuberadius` | `∫₀¹(√(1+a⁴)-1)a⁻²da = 1-2π^{3/2}/Γ(¼)²`; `R=-1-2π^{3/2}/Γ(¼)²`; all 19 printed decimals of `R` | Holomorphy, the pole and its residue | PROVED 4 |
| U02 | `U02_OmegaNonzero` | Prop `prop:an` | `V_n` is an odd integer for every `n`; so `ω_n≠0` for every odd `n≥3`; `ω₁=0` | That `ω_n` is the `ζ(n)` coefficient | PROVED 10 |
| U03 | `U03_R1Positive` | Thm `thm:poly` | `R₁>0` for every `m≥3`; `R₁=0` at `m=2`; `P(1)>0`; so `P(1)R₁=0` only for `m=2` | The identity `P'(1)=P(1)R₁` | PROVED 6 |

Totals: 8 files, 57 statements, all 57 proved, 5231 lines of Lean. `lake build` of the whole library passes
(18296 jobs). No statement of the paper was found false.

So two results that were listed as out of reach are now machine-checked: the closed form of `M₁` (from its
integral `eq:M1t`) and `M₁(α)` at every angle (`thm:Malpha`, `prop:M1short`).

## Wave 3: self-contained analysis theorems

Files are in `lean/Zeta3/Wave3/`. Status from `tools/gate.py check`, run by the planner on 2026-11-25.

| # | File | Paper result | What Lean proved | What it does not cover | Status |
|---|---|---|---|---|---|
| W01 | `W01_GraphLength` | Thm `thm:fn` | The theorem in full, for natural `n`: `|2-L_n-(log n+log f'(1)+1-log 2)/(f'(1)(n-1))| ≤ C log²n/n²` for `n≥N`, with `L_n=∫₀¹√(1+((fⁿ)')²)`; and the five lemmas of the supplement proof (the density `φ`, its moments with `E log S=log 2-1`, the layer-cake formula, the shape and window lemmas) | Real exponents `n`; the higher-order expansion and cumulants after the theorem; that the integral is the arc length | PROVED 6 |
| W02 | `W02_Holo` | Lemma `lem:holo`, eq. `eq:F`, `eq:expansion`, `eq:M1t` | `∫ε(a+b)=2`; `P(p)/4=2-εF(ε)` for `p>1`; the bound `2^{2-Re ε}`; `F` extends holomorphically to `{Re ε<1}`; `P(p)/8-1=Σ_{k≥0}M_k/p^{k+1}` converges for every `p>1` with `M_k=-F^{(k)}(0)/(2k!)`; `M₁` is the integral `eq:M1t` | That the integral `eq:Pint` is the perimeter (W04) | PROVED 7 |
| W03 | `W03_M1Full` | Thm `thm:M1`; the display of `M₀` | `M₁`, the coefficient of `1/p²` in the expansion of `P(p)/8`, equals the closed form of `thm:M1`; `8M₁` is within `10⁻²⁰` of the printed `2.02095276965635835695`; `M₀=λ-√2L` and its printed decimals to `10⁻²⁰`; for every `p>1`, `P(p)/8=1+M₀/p+M₁/p²+Σ_{k≥2}M_k/p^{k+1}` with the two closed forms | `M_k` for `k≥2` | PROVED 6 |
| W04 | `W04_PerimeterIntegral` | eq. `eq:Pint` and the three checks after it | For `p≥1`: the quarter arc is a bijection onto the first-quadrant part of `|x|^p+|y|^p=1`; the four reflected arcs are a simple closed parametrization of the curve; its Euclidean length (total variation) is `(4/p)∫₀¹√(t^{2/p-2}+(1-t)^{2/p-2})dt`; `P(1)=4√2`; `P(2)=2π`; `P(p)→8` | A parametrization-free notion of perimeter | PROVED 11 |

Totals: 4 files, 30 statements, all 30 proved, 3556 lines of Lean. `lake build` of the whole library passes.

Taken together, W04, W02, W03 and D03 give Theorem `thm:M1` from the curve itself: the length of
`|x|^p+|y|^p=1` is `P(p)`, `P(p)/8` has a convergent expansion in `1/p` for `p>1`, and its first two coefficients
are `λ-√2L` and the closed form of `thm:M1`.

## Wave 4: results where ζ(3) itself appears, and three more analysis claims

Files are in `lean/Zeta3/Wave4/`. Status from `tools/gate.py check`, run by the planner on 2026-11-25.
`ζ(3)` and `ζ(5)` are the real parts of Mathlib's `riemannZeta 3`, `riemannZeta 5`.

| # | File | Paper result | What Lean proved | What it does not cover | Status |
|---|---|---|---|---|---|
| X01 | `X01_RhombusDerivs` | Thm `thm:rhombus`; Prop `prop:an` (`a₀…a₃`) | The theorem in full, for all `A,B>0`: the perimeter integral is the true length of `|x/A|^p+|y/B|^p=1` for every `p>0`; `P` is three times continuously differentiable at `p=1`; `P(1)=4r`; the first derivative is `0`, the second `4π²A²B²/(3r³)`, the third `(24A²B²/r³)(ζ(3)-π²/3)`; the log integrals that produce `π²` and `ζ(3)` | The `O(d⁴)` remainder of the boxed expansion | PROVED 17 |
| X02 | `X02_PiP` | Thm `thm:pip` | `π_p=4·2^{-1/p}+E(p)` for every `p>1`; `|π_p-4(1+N₀/p+N₁/p²+N₂/p³+N₃/p⁴)| ≤ C/p⁵` with the printed `N₀…N₃` (so `ζ(3)/4` in `N₃`) | The all-order expansion, its divergence, the growth of `N_k` | PROVED 8 |
| X03 | `X03_RulerConvergence` | Prop `prop:ruler`, last sentence | For every real `q≥1`: `F_q` extends holomorphically to `{Re ε<1}`; the `1/p` series of `L_q(p)/8-1` converges to it for every `p>1` | That the integral `L_q` is the `ℓ^q` length of the curve for `1<q<∞` | PROVED 7 |
| X04 | `X04_PolyFirstDeriv` | Thm `thm:poly` | For every `m≥2`: the boundary of `K_p` is a simple closed curve of finite length `P(p)` for `p` near 1 (both sides); `K₁` is the regular `2m`-gon; the derivative of `P` at `1` is `P(1)R₁`; it is `0` only for `m=2`; it is `2√3 log 2` for the hexagon | Higher derivatives; other polygons | PROVED 15 |
| X05 | `X05_DiamondA4A5` | Prop `prop:an`, Thm `thm:pw` (series through `τ⁵`) | `P` is smooth at `p=1`; `a₀…a₅` as true derivatives, with `a₄=-3ζ(3)-11π⁴/480+π²/2` and `a₅=-7π²ζ(3)/24+6ζ(3)-3ζ(5)/4-2π²/3+11π⁴/120` | `a₆`, `a₇`; the all-order weight claims | PROVED 11 |
| X06 | `X06_Xnyn` | Thm `thm:xnyn` | The theorem in full, for real `n≥2`, with `A(n)` the area integral over the unit square: `|A(n)-(1+√2/n+R(n)/n²+R₃(n)/n³)| ≤ 45(1+log n)⁴/n⁴` with the six printed constants (`d₃` contains `-ζ(3)/2`); the supplement's lemmas (decomposition, Taylor coefficients of `E_g`, Mellin transform, moment identity) | The printed numerical illustrations | PROVED 16 |

Totals: 6 files, 74 statements, all 74 proved, 12668 lines of Lean. `lake build` of the whole library passes
(18315 jobs). No statement of the paper was found false.

Waves 1 to 4: 32 files, 227 statements, 26104 lines of Lean.

## Wave 5: results first listed as OUT

Files are in `lean/Zeta3/Wave5/`. Status from `tools/gate.py check`, run by the planner on 2026-11-25.

| # | File | Paper result | What Lean proved | What it does not cover | Status |
|---|---|---|---|---|---|
| Y01 | `Y01_M2` | Thm `thm:M2` | The theorem in full: `M₂` (the coefficient of `1/p³`, as defined in W02) equals the closed form (T1) of Table `tab:closed`, and `|M₂+0.02216191108726939942|<10⁻²⁰`. The chain: `M₂=-½∫φ₂`; the Euler substitution; `M₂=W₃+T`; `T=(π²+6λ²)/24=(ζ(2)-Li₂(½))/2`; `W₃` equals the raw table value (43 elementary integrals); its reduction with the dilogarithm relations, (I3), (Tr), (D3) and Spence–Kummer. `Li₃` is an integral with its series and `Li₃(1)=ζ(3)` | The sentence that the weight-three part "is a combination of weight three" (it is visible only in the formula); the remarks after the theorem (`Li₃(½)` basis, the 290-digit search, `eq:chi3`) | PROVED (13) |
| Y02 | `Y02_P2` | Thm `thm:P2`, Thm `thm:P2asym` | For `0≤θ<π`: the defect formula, the expansion at infinity, the Mellin identity `∫(√(1+2cs+s²)-s-c)s^{ε-1}ds=(π/sin πε)Λ_ε(θ)`, the exact length formula as a convergent series for `ε=1/(n-1)∉ℕ`, and the removable singularities as limits. For `0<θ<π`: the ranges of `a₁,a₂`, the `ζ(3)` identity for `S_{1,2}`, Taylor bounds that fix `a₁,a₂,a₃`, the uniform bound `|E_n|≤C sin²(θ/2) log n/n³` for `n≥9`, and `thm:P2asym` with error `C log n/n³` | `P_ν` for non-integer `ν` (Mathlib has none; `Λ_ν` is the paper's hypergeometric form); the `n⁻³`, `n⁻⁴` terms; the explicit `9000z/n³` bound; continuation in complex `n` | PROVED (14) |
| Y03 | `Y03_GraphAngle` | Thm `thm:fnangle` | The theorem in full for `0<θ<π`: the density `φ_θ`, its log moment `k₁(θ)=1+cot²(θ/2)log cos²(θ/2)`, the layer-cake identity at angle `θ`, `2-L_n=(1-cos θ)E[m_y(S_θ)]`, and the expansion with error `C log²n/n²`. At `θ=π/2` it is W01 | The Mellin transform of `S_θ`; `k₂`, `k₃`; the higher expansion after the theorem; uniformity in `θ`. Length is the speed integral, as in W01 | PROVED (6) |
| Y04 | `Y04_Fourier` | Thm `thm:fourier` | The theorem in full: for every `k` there are `C_k`, `τ₀` with `|a_k(h_p)-(-1)^k(√2/π)[-4/(16k²-1)+ζ(2)τ²+ζ(3)τ³-(33/16+7k²)ζ(4)τ⁴]| ≤ C_k τ⁵`; with it the diamond coefficients, the layer-integral reduction, the atom integrals, Watson's expansion at `N=4`, and `J₁,J₂,J₃` | `prop:degree`, `prop:theta`; Cauchy's formula `P=∫h_p`; holomorphy in complex `τ` | PROVED (8) |
| Y05 | `Y05_Octahedron` | Thm `thm:oct` | For `S(p)` defined as the paper's face integral `eq:pw:Sn` (`n=3`): integrability for every `p>0`; the graph-area identity for the first-octant piece; the Dirichlet–Gamma factorisation; `E[Y]`, `E[Y²]`; the `d`-series through `d⁴` with error `O(d⁵)`; the log series through `τ⁶` with error `O(τ⁷)`, with `-160ζ(5)/27` at `τ⁵` and `16ζ(3)²/81-572π⁶/76545` at `τ⁶`; both signs | That the face integral is the Hausdorff area (8-fold symmetry, null planes); the weight claims for every `k`; dimensions `N≠3`; `S(2)=4π` and the cube limit | PROVED (11) |
| Y06 | `Y06_RadiusZero` | Thm `thm:zero` | The theorem in full for every `A,B>0`: `limsup|I_n/n!|^{1/n}=1/ρ` with `ρ=√(log²(B/A)+π²/4)`; the Taylor series of `I` and of `P` (in `η` and in `d=p-1`) diverge at every non-zero point; `P` is `C^∞` at `p=1` and not analytic there. With it the Laplace form, the coefficient formula, the zeros of `R₁,R₂`, the Cauchy bound and the exact radius of `H`, the rest bound | The continuation of `H` to the slit disk and the list of zeros in `D(0,ρ₂)` (only their consequences); explicit `C_σ`, `C'_σ`; the text after the theorem (law of `f_n`, Borel, `thm:one`) | PROVED (11) |

Totals: 6 files, 63 statements, all 63 proved, 19294 lines of Lean. `lake build` of the whole library passes
(18328 jobs). No statement of the paper was found false.

Waves 1 to 5: 38 files, 290 statements, 45398 lines of Lean.

## Wave 6: the radius, one more Fourier order, the nonlocal values, the superellipsoid

Files are in `lean/Zeta3/Wave6/`. Status from `tools/gate.py check`, run by the planner on 2026-11-26.

| # | File | Paper result | What Lean proved | What it does not cover | Status |
|---|---|---|---|---|---|
| Z01 | `Z01_RadiusOne` | Thm `thm:one` | The theorem in full: `limsup|M_k|^{1/k}=1`, the radius of `ΣM_kεᵏ` is exactly 1, the sum continues holomorphically across every point of the unit circle except `ε=1`, and not across `ε=1`. With it `P(p)=8(1+M(1/p)/p)` for `p≥1`, and the step from `thm:zero` (no function analytic at `p=1` agrees with `P` on `p≥1`) | Nothing of the theorem. `thm:zero` is imported from Y06 | PROVED (8) |
| Z02 | `Z02_Theta5` | Prop `prop:theta` (`Θ₅`), Prop `prop:degree` (`n≤5`) | `thm:fourier` one order further: `|a_k(h_p)-(-1)^k(√2/π)[…+Θ₅(k²)τ⁵]|≤C_kτ⁶` with `Θ₅(u)=-(¾+4u)ζ(5)-(7/4+4u)ζ(2)ζ(3)`; Watson's expansion at `N=5`; `g_{k,3}`; the weight-five constants `D_{5,0}`, `D_{5,1}` proved analytically; polynomials `Θ₂…Θ₅` of degree `⌊n/2⌋-1` with the stated leading coefficients | `prop:degree` for `n≥6` and the alternating-MZV clause; `Θ₆…Θ₉`; the iterated-integral form of the atoms and the certificates | PROVED (8) |
| Z03 | `Z03_NonlocalValues` | Prop `prop:nonlocal` (last sentence), Lemma `g4:lem:NK` | `ΣH_n(-1)ⁿ/n²=-⅝ζ(3)`, `Σ(-1)ⁿ/n³=-¾ζ(3)`, `ΣH_n2⁻ⁿ/n²=ζ(3)-π²log2/12`, `Σ2⁻ⁿ/n³=⅞ζ(3)-π²log2/12+log³2/6`; Landen's identity for `Li₃`; the Nielsen–Kölbig identity for every `x<1`; `ΣH_nxⁿ/n²=Li₃+S_{1,2}` for `|x|≤1`; the regrouping step of the proof | The formula for `N_e` (a statement about the computer-generated expression of `thm:hexagon`); `Φ_A`, `Φ_L`, the edge data, the hexagon example | PROVED (17) |
| Z04 | `Z04_SuperellipsoidArea` | Thm `thm:S` | For `p≥1`: `S(p)=6G(ε)(P(p)-4)-8ε²Q(ε)` with `Q` as in the theorem and every integral convergent; `S(p)→24`; `|S(p)-24-48M₀/p|≤C/p²`. With it the Gaussian linearisation, Lemma `g3:lem:smooth` at `j=0`, and both parts of Lemma `g3:lem:exact` | That the face integral is the Hausdorff area (as in Y05); the expansion for general `K` and "no `log p` terms"; smoothness of `Q` beyond `j=0`; `s₂`, `s₃`, `Q₀`, `Q₁`; `S_n(p)`; `prop:cuberadius` | PROVED (8) |
| Z08 | `Z08_PrintedDigits` | printed decimals of `M₀`, `κ`, `M₂` | Each printed 20-digit decimal is the truncation of the true coefficient (`W02.Mk`), which is real: `printed ≤ x < printed+10⁻²⁰` for `κ`, `printed-10⁻²⁰ < x ≤ printed` for `M₀`, `M₂` | The other printed decimals of the paper | PROVED (4) |
| Z05 | `Z05_M3` | Thm `thm:M3` | The theorem: `M₃` (the coefficient of `1/p⁴`, as defined in W02) equals `W₄+W₃'+T`, with `W₄` the printed integral, `W₃'` the closed form (T2) and `T=(π²+6λ²)/24`. The chain: `M₃=-½∫φ₃`; the Euler substitution; two integrations by parts; `W₃'` as an integral equals (T2). `Li₄` as an integral with its series. Also: if `W₄=𝓡` then `M₃=𝓡+W₃'+T` | `thm:W4R` itself, `W₄=𝓡`: the paper's proof is a computer certificate over 5184 words, not reproduced; `thm:W4`; the printed decimals of `M₃`, `W₄` | PROVED (11) |
| Z06 | `Z06_P2asymNext` | Thm `thm:P2asym` (next orders) | For `0<θ<π`, real `n≥9`: `|E_n+2z a₂ log n/n³|≤9000z/n³`; the `n⁻³` term with `ζ(3)` in `d(θ)` and error `3.1·10⁴ z log²n/n⁴`; the `n⁻⁴` term with error `1.2·10⁵ z log³n/n⁵`; the supplement's Prop. `g3:prop:n4(a)`. At `θ=π/2`: `a₁`, `a₂`, `a₃=ζ(3)/4-⅓` | The comparisons with Ahmed's paper; the sentence "the logarithm is really there" as a separate statement | PROVED (12) |
| Z07 | `Z07_M2alpha` | Thm `thm:M2alpha`, Thm `thm:flat` (`k=2`) | For `0<α<π` the `u²` coefficient of `M(u;α)` (D04's integral) is `ζ(2)Φ(c)-½Li₂(c²)+Y(c)+Y(-c)`, with `Y` built from harmonic polylogarithms defined as iterated integrals; `Φ(c)=Σc²ⁿ/(4n²-1)` in closed form; `M₂(α)=Σ(r_n+ζ(2)/(4n²-1))c²ⁿ` with the supplement's rational `r_n`; at the square `W₂=T` | The sentence that `W₃` is a combination of weight-three iterated integrals (the 105-integral form); analyticity of `M(u;α)` in `u` (`M₂` is a Peano coefficient); `thm:flat` for `k≠2` | PROVED (10) |

Totals: 8 files, 78 statements, all 78 proved, 13322 lines of Lean. `lake build` of the whole library passes
(18344 jobs). No statement of the paper was found false. In Z05 the planner changed the audited statement block once
after the audit: the unproved `W₄=𝓡` was removed and "Hence `M₃=𝓡+W₃'+T`" became an implication from it; the first ten
statements and all definitions are as audited. The audited original is kept outside the library.

Waves 1 to 6: 46 files, 368 statements, 58720 lines of Lean.

## Wave 7: Cauchy's formula, one more Fourier order, area values, the flat angle

Files are in `lean/Zeta3/Wave7/`. Status from `tools/gate.py check`, run by the planner on 2026-11-26.

| # | File | Paper result | What Lean proved | What it does not cover | Status |
|---|---|---|---|---|---|
| V01 | `V01_CauchyFormula` | Cauchy's formula (Section on Fourier coefficients), Thm `thm:pw` (displayed series) | For `p>1`: `h_p` is the support function of the curve (a maximum over the curve); `P(p)=∫₀^{2π}h_p=π a₀(h_p)` for the true perimeter; the same at the diamond `p=1`. With Z02 at `k=0`: `|P/√2-(4+ζ(2)τ²+ζ(3)τ³-(33/16)ζ(4)τ⁴-(¾ζ(5)+7/4ζ(2)ζ(3))τ⁵)|≤Cτ⁶` for the true perimeter | The all-order and weight claims of `thm:pw`; the rhombus and `N`-dimensional parts; the general convex-curve formula | PROVED (9) |
| V02 | `V02_Theta6` | Prop `prop:theta` (`Θ₆`), Prop `prop:degree` (`n=6`) | `thm:fourier` through `τ⁶` with error `C_kτ⁷`, `Θ₆(u)=(4039/256+733u/8+31u²)ζ(6)-(7/8+2u)ζ(3)²` (the supplement's formula); Watson at `N=6`; `g_{k,4}`; the three weight-six table rows; degree exactly 2 with leading coefficient `31ζ(6)` | `Θ₇…Θ₉`; `prop:degree` for `n≥7` and the alternating-MZV clause; the certificates | PROVED (8) |
| V04 | `V04_AreaValues` | values at `p=1,2` (text before Thm `thm:S`) | `S(1)=4√3` (octahedron) and `S(2)=4π` (sphere) for the face integral of Y05, with the Dirichlet integral `Γ(½)³/Γ(3/2)=2π`. Three more statements (`∫dt/√(t(1-t))=π`, `P(2)=2π`, `P(p)→8`) repeat W04 and are not counted | The enclosed area `A(p)`; `S_n(1)` for `n≠3` | PROVED (4 new of 7) |
| V05 | `V05_FlatAngle` | Thm `thm:flat` (`k≤3`) | For `0<α<π`: `M(u;α)` is smooth at `u=0`; `M_k(α)=Σμ_{k,n}c²ⁿ` for `k=0,1,2,3` with `μ_{0,n}`, `μ_{1,n}` explicit rationals, `μ_{2,n}=r_n+ζ(2)/(4n²-1)`, `μ_{3,n}=r'_n+q'_nζ(2)-ζ(3)/(4n²-1)`, all `r_n,r'_n,q'_n` rational and independent of the angle; `Φ(c)` in closed form; Lemmas `g6:lem:series`, `g6:lem:En`, `g6:lem:gamma`, `g6:prop:formula`, `g6:lem:first` for `k≤3` | `k≥4` and the general weight claim; `Ψ` as a formal series; the uniqueness remark | PROVED (12) |

Totals: 4 files, 33 new statements (36 gated), all proved, 5582 lines of Lean. `lake build` of the whole library
passes (18352 jobs). No statement of the paper was found false. A fifth target, the residue of `prop:cuberadius`,
was proved too but repeats U01 of wave 2 (the planner had not checked); it is kept outside the library and not
counted.

Waves 1 to 7: 50 files, 401 statements (plus 3 repeats), 64302 lines of Lean.

## Wave 8: every order

Files are in `lean/Zeta3/Wave8/`. Status from `tools/gate.py check`, run by the planner on 2026-11-26.

| # | File | Paper result | What Lean proved | What it does not cover | Status |
|---|---|---|---|---|---|
| Q01 | `Q01_PureWeight` | Thm `thm:pw` | To every order, for the squircle and for every rhombus: one sequence `w_k`, each a rational polynomial in `ζ(2),…,ζ(k)` homogeneous of weight `k` (in `ℚ[q²]` for the rhombus), with `|log(P/4r)-Σ_{k=2}^K w_kτᵏ|≤Cτ^{K+1}` for every `K`; the "Hence" for `P/√2` with pure-weight coefficients; `w₂…w₅` explicit; the supplement's lemmas (`lem:pw:pure`, the remainder bound, `lem:pw:asymp`) | The `N`-dimensional surface-area part; the formal-series level (weight is "has a representation of that weight" for the real number); `w_k` for `k≥6` | PROVED (10) |
| Q02 | `Q02_PiPAllOrders` | Thm `thm:pip` | The expansion of `π_p` to every order: `|π_p/4-1-Σ_{k≤K}N_kp^{-k-1}|≤Cp^{-K-2}` for every `K`, with `N_k` by the supplement's formula; `N₀…N₃` agree with X02; the printed closed forms of `N₄` and `N₅`; the Laplace form (`lem:pip:laplace`); `J_{1,j}=j!(1-2^{-j-1})ζ(j+2)`; `lem:pip:C45`; the `h=0` Hankel bound | The growth law of `N_k` and the divergence: stated in Q06, not proved yet | PROVED (13) |
| Q03 | `Q03_Theta7` | Prop `prop:theta` (`Θ₇`), Prop `prop:degree` (`n=7`) | `thm:fourier` through `τ⁷` with error `C_kτ⁸`, `Θ₇` as in the supplement (with `ζ(7)`, `ζ(2)ζ(5)`, `ζ(2)²ζ(3)`); Watson at `N=7`; `g_{k,5}`; the three weight-seven rows; degree exactly 2 | `Θ₈`, `Θ₉`; `prop:degree` for `n≥8` and the alternating-MZV clause | PROVED (9) |
| Q04 | `Q04_FlatAllK` | Thm `thm:flat`, Thm `thm:zladder` (first rung) | `thm:flat` in full: for every `k` and `0<α<π`, `M_k(α)=Σμ_{k,n}c²ⁿ` with `μ_{k,n}` a formal polynomial with rational coefficients in `ζ(2),ζ(3),…` of weight at most `k`, independent of `α`. For every `k`, `n`: the coefficient of `ζ(3)` in `μ_{k+1,n}` is minus the coefficient of `ζ(2)` in `μ_{k,n}`. With it: `log Γ(1+x)+γx=Σ(-1)^wζ(w)x^w/w`, the polygamma values at integers, Lemmas `g6:lem:En`, `g6:lem:gamma`, `g6:prop:formula` for every `k`, and `k=4,5` explicitly | The general rungs of `thm:zladder` (every zeta product, the `ζ(5)` example); `Ψ` as a formal series; the uniqueness remark | PROVED (10) |
| Q05 | `Q05_Decay` | Thm `thm:decay` (lemmas) | Lemma `g5:lem:kernels` (`H`, `K`, their bounds, the residues of `K`); the logit form of `M(ε)` (`g5:lem:logit`); `g5:lem:main`; `g5:lem:F`; `sup u³e^{-(2-√π)u}<114.1`; the arithmetic `<1.21·10⁹`; `√π|ϖ₁|≤1.494`; `C_⋆<10,600,000` | The kernel identity, the representation of `M_k`, `g5:lem:W`, the rest bound, the main-term bound and the theorem itself: stated in Q07, not proved yet. The `254.23` bound of the paper's proof idea (an interval sweep, code only) is not stated | PROVED (10) |

Totals: 5 files, 52 statements, all 52 proved, 11344 lines of Lean. `lake build` of the whole library passes
(18362 jobs). No statement of the paper was found false. Q02 and Q05 came back with some statements unproved; the
planner had the proved statements kept (text unchanged) and the unproved ones moved to Q06 and Q07, which are not
in the library until they are proved.

Waves 1 to 8: 55 files, 453 statements (plus 3 repeats), 75646 lines of Lean.

## OUT: not checkable in Lean with today's Mathlib

These need polylogarithms of weight ≥3, harmonic polylogarithms, multiple zeta values, symbols, Goncharov's
theorem, or certified 300-digit arithmetic. They stay on the written proofs and independent re-checks listed in
the verification map at the top of this repository.

`thm:W4`, `thm:W4R` (`W₄=𝓡`), `thm:wt2`, `thm:pw` (the `N`-dimensional part), `prop:degree` (`n≥8`),
`prop:theta` (`Θ₈`, `Θ₉`), `thm:decay` (in progress, Q07), `thm:varpi`, `thm:pip` (growth of `N_k` and divergence: in progress, Q06),
`thm:Malpha` (`M_k(α)` for `k≥2`; the `M₁(α)` formula is D04), `conj:field`,
`prop:M3angle`, `thm:M4alpha`, `thm:zladder` (the general rungs), `thm:hexagon`, `thm:hexagon5`, `prop:octagon5depth`,
`thm:S` (general `K`, no `log p` terms), `thm:oct` (the weight claims for every order, other dimensions).

This list is what has not been tried or is known to need missing library. Waves 2 to 8 showed that results first
put here can be reached, so an entry here is not a proof that Lean cannot do it.

## Log

- 2026-11-25: checklist written. Lean 4.36.0 and Mathlib 4.36.0 installed; project created in `lean/`; gate
  script `lean/tools/gate.py` written and tested; wave 1 launched (14 targets).
- 2026-11-25: wave 1 back. 13 files pass the gate; T13 has one statement open (second attempt running).
- 2026-11-25: T13 closed by an Opus agent (Bailey's summation without integrals). All 14 files pass the gate;
  `lake build` passes (18280 jobs).
- 2026-11-25: wave 1 statements reviewed by non-Claude models (Codex GPT-6, Kimi K3, GLM-5, DeepSeek on all 14
  files; Gemini on T01 to T10). No statement had to change.
  Thin statements it pointed out: T10 `stmt_4` is true by unfolding (a transcription check of printed constants);
  T12 `stmt_2`, `stmt_4` re-read the display (the content is in `stmt_1`, `stmt_3`).
- 2026-11-25: wave 2 back. All 8 files pass the gate; `lake build` passes (18296 jobs). Wave 3 launched (W01, W02).
- 2026-11-25: wave 3a back. W01 (`thm:fn`, 1506 lines) and W02 (`lem:holo`, 739 lines) pass the gate; `lake build`
  passes (18300 jobs). With D03 this gives `thm:M1` from the perimeter integral. Wave 3b launched (W03, W04).
- 2026-11-25: wave 3b back. W03 (`thm:M1` in full, `M₀`) and W04 (`eq:Pint` is the length of the curve) pass the
  gate; `lake build` passes.
- 2026-11-25: wave 4 back. X01 to X05 pass the gate; `lake build` passes (18314 jobs). X06 (`thm:xnyn`) has four
  statements open; a second run is working on them.
- 2026-11-25: waves 2, 3 and 4 (X01 to X05) reviewed by non-Claude models (Codex GPT-6, Kimi K3, GLM-5, DeepSeek,
  Kilo's free model; Codex ran numerical spot checks). No statement had to change.
- 2026-11-25: X06 (`thm:xnyn`) closed by a second Opus run (3575 lines). All 32 files of waves 1 to 4 pass the
  gate; `lake build` passes. Waves 2 to 4 (without X06) pushed to the public repo in commit `48089b01`. Wave 5
  launched (Y01 to Y06).
- 2026-11-25: wave 5 back. Y02 (`thm:P2`, `thm:P2asym`), Y03 (`thm:fnangle`), Y04 (`thm:fourier`), Y05 (`thm:oct`)
  and Y06 (`thm:zero`) pass the gate; `lake build` passes (18326 jobs). Y01 (`thm:M2`) has 7 of 13 statements open;
  a second run splits them into parts proved in parallel. Statements of X06 and wave 5 sent to
  the non-Claude reviewers.
- 2026-11-25: Y01 (`thm:M2`) closed. Seven Opus agents proved seven parts at the same time in separate files and
  one merged them (6231 lines). The planner's gate passes on all 38 files; `lake build` passes (18328 jobs).
- 2026-11-25: X06 and wave 5 reviewed by non-Claude models (Codex GPT-6 with numerical checks, Kimi K3, GLM-5,
  DeepSeek, Kilo's free model). No statement had to change. Wave 6 launched (Z01 `thm:one`, Z02 `prop:theta` Θ₅, Z03 values of
  `prop:nonlocal`, Z04 `thm:S`).
- 2026-11-25: X06 and wave 5 pushed to the public repo in commit `bc6abbd1` (38 files, 290 statements).
- 2026-11-25: wave 6b launched (Z05 `thm:M3` as `W₄+W₃'+T`, Z06 the next order of `thm:P2asym`, Z07 `thm:M2alpha`).
- 2026-11-25: paper draft 59 (four truncated decimals). No Lean statement changes; line numbers of `main.tex` are
  unchanged.
- 2026-11-25: wave 6c (Z08, one file): the printed decimals of `M₀`, `κ` (draft 59 digits) and `M₂` are
  truncations of the true coefficients. PROVED (4), gate run by the planner, `lake build` passes (18330 jobs).
- 2026-11-26: wave 6 back: Z01 (`thm:one`), Z02 (`Θ₅`, `prop:degree` to `n=5`), Z03 (values of `prop:nonlocal`),
  Z04 (`thm:S`) all PROVED; planner re-ran the gate on each; `lake build` passes (18338 jobs, 43 files, 335
  statements). Non-Claude review of the five finished files started.
- 2026-11-26: wave 7 launched (files in `lean/Zeta3/Wave7/`): V01 Cauchy's formula `P=π a₀(h_p)` and the displayed
  series of `thm:pw` through `τ⁵` for the true perimeter; V02 `Θ₆`; V03 the residue of `prop:cuberadius`; V04
  `S(1)=4√3`, `S(2)=4π`, `P(2)=2π`; V05 `thm:flat` for `k≤2` (stretch).
- 2026-11-26: wave 6b back: Z06 (next orders of `thm:P2asym`) and Z07 (`thm:M2alpha`, `thm:flat` at `k=2`) PROVED.
  Z05 proved `thm:M3` (ten statements) and left `W₄=𝓡` (`thm:W4R`) open: its proof in the paper is a computer
  certificate. The planner removed that statement, kept "Hence" as an implication, re-froze, and gated all three.
  `lake build` passes (18344 jobs): 46 files, 368 statements.
- 2026-11-26: wave 7 back, all five targets proved in about 70 minutes; planner re-ran the gate. V01, V02, V04, V05
  are in the build (18352 jobs): 50 files, 401 statements. V03 repeated U01 and was set aside.
- 2026-11-26: wave 8 launched (files in `lean/Zeta3/Wave8/`, all stretch): Q01 `thm:pw` to every order with the
  weight claim; Q02 `thm:pip` to every order, its divergence and the growth of `N_k`; Q03 `Θ₇`; Q04 `thm:flat` for
  every `k` and the first rung of `thm:zladder`; Q05 `thm:decay`. Wave 7 files sent to the non-Claude reviewers.
- 2026-11-26: waves 6 and 7 reviewed by non-Claude models (Codex GPT-6 with numerical checks, Kimi K3, GLM-5,
  DeepSeek, Kilo's free model, Gemini through `agy`, which works again). No statement had to change. Two comments in
  the header of Z05 were corrected (a draft-58 decimal, a statement number) and the file re-frozen.
- 2026-11-26: waves 6 and 7 pushed to the public repo in commit `70cd9905` (50 files, 404 gated statements, 401
  distinct).
- 2026-11-26: wave 8 back. Q01 (`thm:pw` to every order), Q03 (`Θ₇`), Q04 (`thm:flat` in full, first rung of
  `thm:zladder`) PROVED. Q02 (`thm:pip`) proved 13 of 20 statements and Q05 (`thm:decay`) 10 of 18; each was split
  into a finished file and a continuation (Q06, Q07). Planner gated all five; `lake build` passes (18362 jobs): 55
  files, 453 statements. Continuation run launched for Q06 and Q07 (parts in
  parallel, then merge). Non-Claude review of wave 8 started.
- 2026-11-26: non-Claude review of wave 8 back (Codex with numerics, Kimi, GLM, DeepSeek, Kilo's free model; `agy`
  is out of quota again). No statement had to change. Codex and Kimi each derived `e_{0,5} = 721/3840`, which the
  supplement does not print. One process point: `V05.HrQ` and `V05.rhoQ`, used in the frozen statements of Q04,
  are defined below the END STATEMENTS line of V05; the planner read them, and Lean proves their casts equal the
  audited `V05.Hr` and `V05.rho`. From now on a definition used in a frozen statement must sit in a frozen
  block.
