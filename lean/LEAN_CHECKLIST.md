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
| X06 | `X06_Xnyn` | Thm `thm:xnyn` | in progress | | not finished |

Totals so far: 5 files, 58 statements proved, 9093 lines of Lean. `lake build` of the whole library passes
(18314 jobs). No statement of the paper was found false.

All waves: 31 files, 211 statements, 22529 lines of Lean (X06 not counted).

## OUT: not checkable in Lean with today's Mathlib

These need polylogarithms of weight ≥3, harmonic polylogarithms, multiple zeta values, symbols, Goncharov's
theorem, or certified 300-digit arithmetic. They stay on the written proofs and independent re-checks listed in
the verification map at the top of this repository.

`thm:M2`, `thm:M3`, `thm:W4`, `thm:W4R`, `thm:wt2`, `thm:pw` (all orders), `thm:fourier`, `prop:degree`,
`prop:theta`, `thm:zero` (main claim), `thm:one`, `thm:decay`, `thm:varpi`, `thm:pip` (all orders and growth of `N_k`),
`thm:P2` (main claim), `thm:P2asym`, `thm:Malpha` (`M_k(α)` for `k≥2`; the `M₁(α)` formula is D04), `conj:field`, `thm:flat` (main claim), `thm:M2alpha`,
`prop:M3angle`, `thm:M4alpha`, `thm:zladder`, `thm:hexagon`, `thm:hexagon5`, `prop:octagon5depth`,
`thm:fnangle`, `thm:S`, `thm:oct` (main claim).

This list is what has not been tried or is known to need missing library. Waves 2 to 4 showed that results first
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
