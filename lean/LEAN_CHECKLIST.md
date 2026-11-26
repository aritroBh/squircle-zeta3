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
- Wave 2 builds a small dilogarithm library and checks the dilogarithm identities the M₁ proof rests on.
- The main theorems (M₁ to M₄, pure weight, radius zero, decay law, hexagon theorems) stay verified the way
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

If Lean finds that a statement in the paper is false as written, that is reported to the authors first.

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

## Wave 2: needs a small dilogarithm library first

| # | Paper result | What Lean checks | Status |
|---|---|---|---|
| 19 | library | `Li₂` on `(-∞,1]` by series and integral; reflection, inversion, duplication, Landen | TODO |
| 20 | Thm `thm:M2` eq. (T) | `Li₂(½)=π²/12-½log²2`, so `(π²+6λ²)/24=(ζ(2)-Li₂(½))/2` | TODO |
| 21 | Lemma `lem:dilog` | the dilogarithm identities at `b=√2-1` | TODO |
| 22 | Prop `prop:M1short` vs Thm `thm:Malpha` | at `α=π/2` the four-dilogarithm form equals the `M₁` of Thm `thm:M1` | TODO |
| 23 | Prop `prop:cuberadius` | `∫₀¹(√(1+a⁴)-1)a⁻²da = 1-2π^{3/2}/Γ(¼)²` | TODO |
| 24 | Prop `prop:an` | `ω_n≠0` for every odd `n≥3` | TODO |
| 25 | Thm `thm:poly` | `R₁>0` for `m≥3` (so `P'(1)=0` only for `m=2`) | TODO |

## Wave 3: self-contained analysis theorems (long; only if waves 1 and 2 go well)

| # | Paper result | Note | Status |
|---|---|---|---|
| 26 | Thm `thm:fn` | length of the graph of `fⁿ`; real analysis only, no special functions | TODO |
| 27 | Lemma `lem:holo` | holomorphy of `F` and convergence of the `1/p` expansion | TODO |
| 28 | Thm `thm:M1` | full closed form of `M₁`; needs 19 to 22 and the integral reduction | TODO |

## OUT: not checkable in Lean with today's Mathlib

These need polylogarithms of weight ≥3, harmonic polylogarithms, multiple zeta values, symbols, Goncharov's
theorem, or certified 300-digit arithmetic. They stay on the written proofs and independent re-checks listed in
`VERIFICATION_MAP_draft58.md` at the top of this repository.

`thm:M2`, `thm:M3`, `thm:W4`, `thm:W4R`, `thm:wt2`, `thm:pw` (main claim), `thm:fourier`, `prop:degree`,
`prop:theta`, `thm:zero` (main claim), `thm:one`, `thm:decay`, `thm:varpi`, `thm:pip` (asymptotics),
`thm:P2` (main claim), `thm:P2asym`, `thm:Malpha`, `conj:field`, `thm:flat` (main claim), `thm:M2alpha`,
`prop:M3angle`, `thm:M4alpha`, `thm:zladder`, `thm:hexagon`, `thm:hexagon5`, `prop:octagon5depth`,
`thm:fnangle`, `thm:xnyn` (main claim), `thm:S`, `thm:oct` (main claim).

## Log

- 2026-11-25: checklist written. Lean 4.36.0 and Mathlib 4.36.0 installed; project created in `lean/`; gate
  script `lean/tools/gate.py` written and tested; wave 1 launched (14 targets).
- 2026-11-25: wave 1 back. 13 files pass the gate; T13 has one statement open (second attempt running).
- 2026-11-25: T13 closed by an Opus agent (Bailey's summation without integrals). All 14 files pass the gate;
  `lake build` passes (18280 jobs).
