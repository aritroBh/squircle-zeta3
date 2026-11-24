# Referee report — g8 pack (Theorem 20, Proposition 5, Theorem 21 of the corner paper)

Scope: the three new sections `sec_g8_M2alpha.tex`, `sec_g8_M3alpha.tex`, `sec_g8_M4alpha.tex`, their machine steps M1–M9 and checks, against the paper text (`main_tex_statements_l1500-1580.tex`) and the refereed g6 sections they cite. I read every file in the pack (all .tex, all macros, REPORT_G8.md, wrap_g8.pdf via pdftotext, all scripts and .out files), re-ran every supplied script from a copy in my working directory (all reproduce their .out files exactly), checked the new proofs on paper, and verified by hand the specific items listed below. I did **not** get to write my own independent code for some of the tasked recomputations; those are listed explicitly under "Not checked".

---

## Verdicts

| Claim | Verdict |
|---|---|
| Thm 20, all-n identity (M1/M2 + criterion, `sec_g8_M2alpha.tex:37-124`) | **PROVED** (argument audited; machine step reproduced; only the soundness direction of the normal form is used, and it is valid) |
| Thm 20, classical form (M3 + M4, `:126-224`) | **PROVED** (method sound; script reproduced 33/33; my independent hand-recomputation of 3 of the 10 identities NOT done — see "Not checked") |
| Prop 5 (symbol, δ₂,₂, one-term impossibility, two-term identity, `sec_g8_M3alpha.tex`) | **PROVED AFTER REPAIRS** (mathematics correct and exhaustive over the stated classes; the *paper's* wording needs two fixes — findings 4, 5 — which the supplement itself already demands) |
| Thm 21 (ζ(3) ladder, π⁴ part, M8/M9, counts, `sec_g8_M4alpha.tex`) | **PROVED AFTER REPAIRS** (proof correct and reproduced; the *paper's* "22 such terms of weights one to three" is inaccurate — finding 6) |
| Gap boxes / coverage honesty of the supplement | **Accurate.** The supplement's "what is not proved" boxes match my own assessment of what is and is not proved. The paper overclaims in two places; the supplement discloses both. |

---

## Priority 1 — the all-n argument (`sec_g8_M2alpha.tex:41-121`)

I checked the four lemmas on paper and audited the implementation line by line.

* **Lemma g8:lem:formal (equal elements ⇒ equal values).** Correct: an element of 𝔇 is a finite formal sum Σ R_w[w] with R_w ∈ ℚ(N) whose poles lie at integers ≤ 1; evaluation at N = 2n ≥ 2 never hits a pole. Only this *soundness* direction is ever used by the machine steps (run_k requires the difference dict to be literally empty). The converse (could two different formal elements agree as functions of n?) is **not needed** for the proof; if it failed it would cause false alarms, not false proofs. sympy's `field('n', QQ)` gives canonical rational-function normal forms, so "difference reduces to 0" *does* decide equality of the coefficients.
* **Lemma g8:lem:shift (j ≥ -2).** Correct. S_j(v) recurses on the first letter with denominators (N+t), t ≥ j+1 ≥ -1; the recursion stays in range because S_{t-1} has t-1 ≥ -2. Poles lie at N ≤ 1, hence never at N = 2n ≥ 2. The code `_shift` matches the lemma (including `sa**(t%2)` = σ^t for σ=±1, valid also for negative t).
* **Lemma g8:lem:hcoef (δ ≥ -1).** Correct: [x^{2n+δ}]H_w = cf_w ε^δ (N+δ)^{-κ} ev_n(S_{δ-1}(v_w)) with δ-1 ≥ -2; `hpl_coeff_terms` implements exactly the recursion of Lemma g6:lem:dict (single term per word; cf, ε, κ, v all match my hand derivation), and `hpl_coeff_at` applies the shift.
* **Lemma g8:prop:crit.** Correct and complete: even coefficients N=2n ≥ 2 by (i) L_0 = T_z; odd N=2n+1 ≥ 3 by L_1 = 0; the remaining N ≤ 1 by (ii), an exact finite computation. In the script the low-order check covers N = -4..1, which contains every N that can occur (e ≥ -3 and H_w is supported on N ≥ 1, so c^e H_w is supported on N ≥ -2). It uses the *independent* direct recursion `hpl_coeffs` (verified against the HPL definition by hand), not the shift machinery — so (ii) does not rely on the code being tested in (i). Constant terms (e.g. the 1/2 of Ψ and Ψ₂) are correctly excluded from L_0 and included in the N ≤ 1 check.
* **Denominator audit (`:93`).** From the code: denominators are n, n+t (t ≥ 0, from χ⁽ʳ⁾₁), 2n, 2n-1, 2n+t (t ≥ -1 from shifts/hcoef; t ≥ 0 from ρ_j and χ⁽ʳ⁾₃). Poles at n ∈ {0, 1/2, -1/2, -1, ...} — never at an integer n ≥ 1. The text's parenthetical attribution ("n, n+t and 2n+t with integers t ≥ -1") is garbled (n+t occurs only with t ≥ 0), but the conclusion is right. **Finding 2 (minor).**
* **M1/M2 vs the text.** `m2_terms_from_statement()` encodes exactly the eleven terms (e,w;q) listed at `:107-108`, the even sign s = (-1)^{e+ν(w)} matches Lemma g6:lem:ops(a), and the extra term -1/(2n²)[∅] is exactly -½ Li₂(c²). The check is precisely the criterion's (i)+(ii) for L = Y(c)+Y(-c) - Σ c^{2n}/(2n²). Reproduced: "THEOREM 20 identity (all n): True", "[ζ(2)]P_{2,n} == 1/(4n²-1): True". The mu-side construction matches Definition g6:def:formal + Lemma g6:lem:A (E-series, exp(Δ_⋆), χ-substitution, Λ_j with the t=0 factor 1/2 folded into -1/(n(2n-1))·Π(n+t)/(2n+t) — verified algebraically). The selftest (98 words × n=1..7 × δ=-1..4 = 4116 checks, plus harmonic numbers and stuffle against direct computation) reproduces True.
* **Machine step reproduction** (run in `ref_scripts/`, copies of the pack's scripts): `g8_allN.py selftest|2|3|4` → outputs identical to the supplied .out files, all checks True.

## Priority 2 — classical form (M3, M4)

* Method (verified on paper): prove C_w' = f_{a₁}C_{w'} exactly (chain rule; polylog(1,·) rewritten by hand into {log x, log(1-x), log(1+x), log 2}, all factors positive on (0,1) — I checked all eight entries of `NEG_LOG_1M`); then vanishing is `sp.cancel()==0` on rational-function coefficients of a fixed monomial basis — sufficient, uses no relation among basis elements. Base point: exact limit at x→0⁺ with the proved special values (Li₂(1), Li₂(-1), Li₂(1/2), Li₃(1), Li₃(-1), Li₃(1/2)) substituted before comparison — the "lim C_w = 0" check is exact. Branches: all eight arguments are strictly monotone with ranges in (-1,1] (verified against the printed table at `:130-135`).
* Reproduced `g8_cor3_verify.py` (17.8 s): 33/33 PROVED (6 weight-two identities, 5 dilog reductions + duplication, 10 weight-three identities "derivative identity EXACT True; constant EXACT 0 True", range grid check, sign-rule numeric check, and the 31-coefficient comparison in ℚ(c) all `diff = 0`, ζ(3) coefficient 7/4+2c).
* The corollary's displayed formula (`:190-198`) matches the script's typed list term by term (checked all 31 by reading) — Li₃ at the eight arguments c, -c, 1-c, 1/(1+c), (1+c)/2, ±(1-c)/(1+c), 2c/(1+c), matching the paper's sentence.
* **Observation (finding 8, cosmetic):** the .out prints "monomials checked: 1" for every weight-3 word; the pass/fail flags come from `cancel()==0` over the basis coefficients, and I verified the method is sound, but I did not instrument the code to reconcile that printed count.

## Priority 3 — Proposition 5

* Definitions (V, symbol with outermost letter last, B₂ = T²/shuffles ≅ Λ²V, δ₂,₂ = (x₁∧x₂)∧(x₃∧x₄)) are precise and standard; I verified the differentials of Li_{1,1}, Li_{1,2}, Li_{2,1}, Li_{2,2} by hand (all four correct), the specialisation maps at x=1/2 (correct, including [−c]=[c] and [1∓c)/2]=[1∓c]−[2]), and that with X,U↦0 only (Y∧V₁)∧(Y∧W) survives, giving δ₂,₂ = -p∧r (y=1-c) and +q∧r (y=1+c) — hence the two-term identity with coefficients (1-c)/c, -(1+c)/c equals δ₂,₂(Ω). 
* Lemma g8:lem:vanish (δ₂,₂ kills shuffle products and Li₄): proofs checked by hand, correct; the exhaustive machine controls (900 random + all 243 shuffle products, 16 Li₄ tensors, positive control) reproduce.
* Lemma g8:prop:one (one term cannot): correct and exhaustive over the stated class (r(c)·G, G any rational-constant-coefficient symbol in any letters, constants projected out): δ₂,₂(r·π^{⊗4}τ) has component ratio b/d ∈ ℚ, but Ω's ratio is (1-c)/(1+c) ∉ ℚ. This is exactly the class the paper's sentence names.
* Reproduced `g8_M3_delta22.py`: δ₂,₂(Ω) = (1-1/c)([c]∧[1-c])∧([1-c]∧[1+c]) - (1+1/c)([c]∧[1+c])∧([1-c]∧[1+c]), "equals the displayed formula of Proposition 5: True"; two-term identity residual 0 modulo the constant letter [2]; with [2] kept as a formal letter, 4 nonzero residual components (printed).
* **Finding 4 (must fix, paper wording; already the writer's disagreement #4):** the paper's sentence "Two terms are enough: δ₂,₂(Ω) = δ₂,₂((1-c)/c Li_{2,2}(1/2,1-c) - (1+c)/c Li_{2,2}(1/2,1+c))" is literally true only *modulo the constant letter [2]* (M7 prints 4 residual components with [2] kept). The proposition's own "(of the function, modulo constants)" arguably carries over, but the sentence should say it. Repair: append "(symbol level, modulo constants)" — exactly what the writer recommends; the supplement's gap box (item 5) already discloses the residual.
* **Finding 5 (must fix, paper overclaim):** "So M₃(α) is not classical as a function of the angle" is **not proved** (needs invariance of the symbol under the choice of expression, i.e. independence of HPLs and classical functions over ℚ(c)). The supplement's gap box (item 2) discloses this. Truthful replacement: "So the weight-four symbol of the displayed expression of M₃(α) is not that of a classical combination; concluding that the function M₃(α) is not classical needs the linear independence of its harmonic polylogarithms and classical functions over ℚ(c), which is not proved."
* Note: the claim "checked by two independent programs" is defended by the writer as COR8B + COR8BR (+ their own M7); COR8R could *not* confirm (convention mismatch). I can only repeat the writer's account; my own reproduction of M7 is a third implementation consistent with the identity modulo [2].

## Priority 4 — Theorem 21

* ζ(3) part: [ζ(3)]P_{4,n} = -[ζ(2)]P_{3,n} = -q'_n via Lemma g6:thm:ladderprecise with Q_z = σ₁σ₂, γ₁=1 — verified by hand; M9 confirms the ladder rung as an identity in 𝔇 (reproduced True).
* π⁴ part: partial fractions verified by hand: (n+3)/((2n-1)(2n+1)(2n+3)) = 7/16·1/(2n-1) - 5/8·1/(2n+1) + 3/16·1/(2n+3) and -n/(...) = -1/16·1/(2n-1) - 1/8·1/(2n+1) + 3/16·1/(2n+3) (checked at the three poles); the three sums Σc^{2n}/(2n±1), Σc^{2n}/(2n+3) standard; F₁ constant 5/8-1/16 = 9/16, F₂ constant 1/8-1/16 = 1/16; 32F₁+80F₂ = (9c-30/c+21/c³)artanh c + 23 - 21/c² (all five coefficients checked). `g8_pi4.py` reproduces all True.
* **Finding 6 (must fix, paper wording; writer's disagreement #3):** the paper says "Ψ₂ is a rational combination of 22 such terms of weights one to three". The file has **21 HPL terms (3, 6, 12 of weights 1, 2, 3) plus the constant 1/2**; "22 such terms of weights one to three" is inaccurate. R₄ = 109 = 1/8/27/73 (weights 2/3/4/5) confirmed by M8 output. The supplement states the correct count (`sec_g8_M4alpha.tex:22`); the paper should read "... of 21 terms c^e H_w(±c) of weights one to three, plus the constant 1/2".
* M8 reproduced: parts 1 (109 terms) and ζ(2) (22 terms), even/odd identities and low orders True; [ζ(4)]P_{4,n} = (n+3)/((2n-1)(2n+1)(2n+3)), [ζ(2)²]P_{4,n} = -n/((2n-1)(2n+1)(2n+3)).

## Priority 5 — independent numerics

Reproduced `g8_numeric.py` (6.5 s) exactly: M₂ from the defining Cauchy integral vs Thm 20 (HPL form) and vs the classical form at α = π/3, π/2, 2π/3, 1.8: differences **4.45e-39, 3.38e-39, 2.33e-39, 2.28e-39** (HPL) and **4.31e-39, 3.53e-39, 2.42e-39, 2.39e-39** (classical); M₄ vs Thm 21: **1.83e-37, 1.7e-37, 1.41e-37, 1.5e-37** (two obtuse angles included: 2π/3 and 1.8 rad). I read the script: the Cauchy integral (r=0.2, 64 nodes, 40 digits) and the HPL Taylor evaluations are correct as written; the M₄ check uses the json term lists (the actual data of the proof). **Not done by me:** an independent (own-code) numerical check; **M₃ is not numerically checked anywhere in the pack** (`sec_g8_M3alpha.tex:238` says so) and I did not get to run my own M₃ check — flagged as a gap in my review, not in the supplement's logic (M5 is the proof). Planted-error tests: all 6 mutants in `g8_mutants.py` are caught (reproduced); I did not plant my own.

## Numbered findings

1. **(minor)** `macros_all.tex` / `macros_g6.tex`: `\gapnote` prints "Step not written out in our notes:", but the g8 sections use it for "What the section does not prove" boxes (the REPORT calls them "gap boxes"). Cosmetic but confusing; give the box its own heading.
2. **(minor)** `sec_g8_M2alpha.tex:93`: the denominator attribution is garbled (the n+t family occurs only with t ≥ 0); the conclusion (no pole at an integer n ≥ 1) is correct as verified from the code. Repair the parenthetical.
3. **(note)** Only the soundness direction of the normal form is used (equal formal elements ⇒ equal values for every n ≥ 1); valid. Completeness is irrelevant to the proof. 
4. **(must fix — paper's Prop 5 sentence)** add "modulo the constant letter [2]" to the two-term identity; with [2] kept, 4 residual components remain (verified).
5. **(must fix — paper overclaim)** "M₃(α) is not classical as a function of the angle": only the symbol of the displayed expression is proved non-classical; invariance is not proved (supplement's gap box item 2 says so). Qualify the sentence.
6. **(must fix — paper's Thm 21 wording)** "22 such terms of weights one to three" → 21 HPL terms + the constant 1/2 (supplement states this correctly).
7. **(note)** The mu-side symbolic-n implementation is cross-checked only against the writer's own second code (`g6_common.py`, n=1..8); Lemma g6:lem:A it implements is in the refereed g6 sections. I verified the implementation matches the construction by reading; no independent re-implementation by me.
8. **(note)** `g8_cor3_verify.out`: "monomials checked: 1" printed for every weight-3 word — cosmetic; the PROVED flags derive from exact `cancel()==0` checks whose method I verified.
9. **(note)** External referee reports (COR3k, COR8R/B/BR, G1/G2R) are taken from the writer's summary; I could not audit them.

## What I recomputed / verified (evidence)

- Re-ran in `ref_scripts/` (exact copies): `g8_allN.py selftest` (4116 checks True), `g8_allN.py 2|3|4` (all identities True; counts 11 / 74+15 / 109+22 terms; ladder True; [ζ(4)], [ζ(2)²] formulas as printed), `g8_M3_delta22.py` (δ₂,₂(Ω) = displayed formula True; controls; two-term residual 0 mod [2]; 4 residuals with [2] kept), `g8_cor3_verify.py` (33/33 PROVED, 17.8 s), `g8_pi4.py` (all True), `g8_mutants.py` (6/6 mutants caught), `g8_numeric.py` (residuals listed above). All outputs identical to the supplied .out files.
- By hand: the four lemmas of the criterion; the shift/hcoef recursions against `g8_zsum.py`; the mu-side construction against Definition g6:def:formal/Lemma g6:lem:A (including the Λ_j/ρ_j factor check); the eleven M₂ terms against the statement; the four Li_{a,b} differentials; the specialisation maps; the one-term impossibility ratio argument; the vanishing lemma; the π⁴ partial fractions and the 1/2880 combination; the low-order logic and coverage of all N; the corollary's 31 atoms against the script's typed list; the 8 log-splittings in `NEG_LOG_1M`.

## Not checked

- Independent (own-code) recomputation of ≥3 of the ten weight-three identities (tasked; I verified the method and the script's reproduction only).
- Independent own-code recomputation of δ₂,₂(Ω) from the 56-term data (I read the implementation — it is correct as written — and re-ran it).
- Own numerical evaluation of M₂/M₃/M₄ (I relied on the pack's `g8_numeric.py`, which I read and reproduced; M₃ has no numerical check anywhere).
- Planting my own mutant (the pack's 6 mutants are all caught).
- The g6 sections are taken as already refereed (per instructions); I only verified the cited lemmas' statements and the g8 usage.
- The wrap build (tectonic) was not re-run; I inspected `wrap_g8.pdf` via pdftotext (box headings, no "??").

Bottom line: no fatal errors. The three new sections are mathematically sound and honest about their gaps; Theorem 20 (both forms), the all-n machinery, and Theorem 21's proof check out, and Proposition 5's mathematics is correct — subject to the two wording repairs in the paper's Proposition 5 paragraph (items 4, 5) and the "22 terms" count (item 6), all of which the writer's own report already flags.