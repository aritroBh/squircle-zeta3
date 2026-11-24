# COR4: M₃(α) is classical only at special angles — the generic (2,2) obstruction is nonzero; the square corner is special because Λ²B₂ collapses there

## 1. The question and the answer

F(c) := weight‑4 part of R₃(c) (56 of the 74 terms of `hpl_fit_k3_(0,0,0).json`, weights 4, coefficients in ℚ[c,1/c]).

**Answer (symbol level, PROVED).** F is **not** classical as a function of c. Its antisymmetrised (2,2) component (Goncharov's cobracket δ₂,₂: B₄ → Λ²B₂, symbol level) is nonzero as a function of c:

> δ₂,₂(F) = (1 − c⁻¹)·([c]∧[1−c])∧([1−c]∧[1+c]) + (−1 − c⁻¹)·([c]∧[1+c])∧([1−c]∧[1+c]) ≠ 0

The (3,1) and (1,3) components are also nonzero. Specialising:

| angle α | c | δ₂,₂(F) |
|---|---|---|
| π/2 | 1/√2 | **0** ✓ (required; paper Thm thm:W4R) |
| 2π/3 | 1/2 | 0 (automatic: 2 letters ⇒ Λ²B₂ = 0) |
| (rational) | 1/3 | 0 (automatic: 2 primes) |
| π/3 | √3/2 | **≠ 0**: −2·([2]∧[2+√3])∧([3]∧[2+√3]) + (4√3/3)·([2]∧[3])∧([2]∧[2+√3]) |
| (rational) | 2/5, 3/7, 4/7, 5/9, 5/11 | **≠ 0** (12 components each) |
| 3π/4 | cos(3π/8) | **≠ 0** (3 components, quartic field) |
| π/6 | cos(π/12) | **≠ 0** (3 components, quartic field) |

So M₃(α) is classical only at special angles. The mechanism: δ₂,₂(F) vanishes exactly where the letters collapse so that B₂ becomes 1‑dimensional (then Λ²B₂ = 0). At c = 1/√2 the letters live in span{[2], [1+√2]}, hence Λ²B₂ = 0 automatically — the (2,2) obstruction that is present for generic c disappears at the square corner. (Note: δ₂,₂ = 0 at c = 1/2, 1/3 is likewise automatic and is *not* evidence for classicality there; the finer Li₄‑span test at those points was not completed — see §6.)

Caveat: this is a symbol‑level statement. Vanishing δ₂,₂ is necessary for classicality; non‑vanishing certifies non‑classicality. As for M4c, translating symbol identities into identities of actual numbers is modulo the standard period conjecture (kernel of the full symbol). Label: **PROVED** for the exact symbol algebra, **NUMERICAL** for number evaluations.

## 2. Task 1 — symbol code and its validation (PASS)

Module `cor4_sym.py`. Conventions: H_{a,w}(x) = ∫₀ˣ f_a H_w with f₀=1/t, f₁=1/(1−t), f₋₁=1/(1+t) (as in result.md); S(H_{a₁…a_n}(x)) = g(a_n)⊗…⊗g(a₁) with g(0)=[x], g(1)=−[1−x], g(−1)=[1+x] (verified against the standard S(Li_n(x)) = −[1−x][x]ⁿ⁻¹ via the coproduct Δ̃Li₂ = (−log(1−x))⊗log x). Rational‑function coefficients c^e are carried as scalars (at rational c₀ this is exactly the number‑level symbol; at quadratic points coefficients in ℚ(√d) are carried as (p,q) pairs, as in M4c).

Validation (`step10_symbol_checks.py`, all PASS, exact):
- S(H_{0,0,0,1}) = −[1−x][x]³ = S(Li₄(x)); S(H_{0,0,0,−1}) = +[1+x][x]³ = −S(Li₄(−x)) (note H_{0,0,0,−1} = −Li₄(−x), matching hpl_eval.py's H_{0,0,−1} = −Li₃(−x)).
- Shuffle consistency: S(H_{0,1})² = S(2H_{0,1,0,1} + 4H_{0,0,1,1}).
- S(Li₄(x)+Li₄(−x)) = (1/8)S(Li₄(x²)) using [1−x²]=[1−x]+[1+x].
- δ₂,₂(Li₄)=0; δ₂,₂(Li₃,₁(x,1))=δ₂,₂(Li₂,₂(x,1))=0 (one‑variable depth‑2 is (2,2)‑invisible); positive control δ₂,₂(H_{0,1,0,−1}(−x)) = ((0,1),(0,2)) ≠ 0.
- Numeric (dps 60): H_{0,0,0,±1}(x)=±Li₄(±x), H_{0,0,1,1}=Li₃,₁(x,1), H_{0,1,0,1}=Li₂,₂(x,1), H_{0,0,1,−1}=−Li₃,₁(x,−1), H_{0,0,−1,1}=−Li₃,₁(−x,−1), H_{0,−1,0,1}=−Li₂,₂(−x,−1) at x=0.31, 0.77 — max residual 1.4·10⁻⁶⁰ (my Li_{m,n}(x,y) := Σ_{n>m≥1} xⁿyᵐ/(nᵐ·mⁿ‑pattern): Li₃,₁(x,y)=Σ xⁿyᵐ/(n³m), Li₂,₂=Σ xⁿyᵐ/(n²m²)).

## 3. Task 2 — symbol of F and the (2,2) component (PROVED)

F's symbol as a function of c (56 pure tensors over {[c],[1−c],[1+c]}, coefficients c^e, e∈{−1,0,1}) was computed exactly, printed in full (`step11_F_symbol.py` output), and saved to `F_symbol.json`. δ₂,₂(F) as stated in §1: two nonzero wedge components with Laurent‑polynomial coefficients −c⁻¹+1 and −c⁻¹−1. The (3,1)/(1,3) components are nonzero at each c‑degree (36/26, 42/14, 58/36 terms for e = −1, 0, 1).

**Important correction made en route.** A first version used only (a∧b)⊗(c∧d) (tensor square). That gave a nonzero value at c = 1/√2, contradicting the paper's *proved* theorem (W₄ classical, exact certificate). I verified the pipeline against the paper numerically (below) and traced the discrepancy to the missing second antisymmetrisation: the true Goncharov component targets Λ²B₂ (exterior square), which also kills symmetric terms like (P∧Q)⊗(P∧Q) coming from Li₂(b)². After the fix, the square corner gives exactly 0, as required — a strong end‑to‑end control of the whole pipeline.

Cross‑checks against the paper (NUMERICAL, dps 60–80):
- F(1/√2) = −0.07245446822891105817780684456332961285746… (my fast exact‑series HPL evaluator `hpl_series.py`, cross‑checked vs polylogs).
- Paper's W₄ formula evaluates to 0.795785059832240811320359819629726762546527414… vs the paper's printed 0.795785059832240811320359819629726762546527413871… ✓.
- D = F − W₄ = −0.86823952806115186949816666419305638… equals the explicit decomposable combination ζ(3)Φ − ζ(2)Ψ − R₃^{≤3} + W₃' + T with residual 3.61·10⁻⁵⁹ ⇒ [F(1/√2)] = [W₄] in B₄ (numbers).
- Symbol of the paper's Li₄ part: δ₂,₂ = 0 ✓; symbol of the paper's product terms: δ₂,₂ = 0 ✓ (after the fix).
- M₃(G1 formula) at c₀ = 1/√2 = 0.02277052149178657745673050974232780025018… matches the paper's M₃ (|diff| = 2.61·10⁻⁵²).

## 4. Task 3 — specialisations (table in §1; PROVED at symbol level)

Method: per point, letters are mapped to a verified multiplicative basis of the field (primes for ℚ; {2, 1+√2} for ℚ(√2); {2, 3, 2+√3} for ℚ(√3); {d, 2−d, 2+d} with 2 = d²(2−d)(2+d) for cos(3π/8); {2, d, 2+d} with (2+d)(2−d) = d⁻² for cos(π/12)); coefficients c₀^e carried exactly (pairs / power‑basis tuples); δ₂,₂ computed per c‑exponent with rational arithmetic and combined exactly.
- Field relations verified **exactly** with sympy: 2 = d²(2−d)(2+d) at cos(3π/8) (minpoly x⁴−4x²+2); (2+d)(2−d) = d⁻² at cos(π/12) (x⁴−4x²+1).
- Multiplicative independence of the quartic bases: log‑embedding determinant ≈ −1.973144 ≠ 0 (cos(3π/8)) and rank 3 (cos(π/12)) — NUMERICAL.
- The √3/2 value was reproduced from the function‑level formula by hand: −(A+B) = 4√3/3 on ([2]∧[3])∧([2]∧[2+√3]), −A+B = −2 on ([2]∧[2+√3])∧([3]∧[2+√3]), A = 1−2√3/3, B = −1−2√3/3 ✓ consistent with the printed machine result.
- Results saved to `delta22_specializations.json`.

Interpretation: δ₂,₂(F) is a fixed nonzero tensor; at a specialisation it can only vanish if the letter images collapse, and at every tested point with ≥ 3 independent letters (√3/2, 2/5, 3/7, 4/7, 5/9, 5/11, cos(3π/8), cos(π/12)) it is nonzero ⇒ **M₃ is not classical at π/3, 3π/4, π/6, or generic rational c** (symbol level). The square corner gives 0 because B₂ collapses to one dimension.

## 5. Task 4 — depth‑two part: status (incomplete, stated honestly)

Characterisation established: the δ₂,₂‑visible depth‑two content of F must come from depth‑2 MPLs with *mixed* letters (one‑variable Li₃,₁(x,1) and Li₂,₂(x,1) have δ₂,₂ = 0; F's obstruction comes from pure tensors like H_{0,1,0,−1}(c) = −Li₂,₂(c,−1) whose letters straddle [1−c] and [1+c]). Built for the decomposition: 40 S‑unit arguments r(c) with r, 1−r factoring over {2, c, 1−c, 1+c} (`step13_library.py`, `sunit_args.json`), the candidate dictionary H_{0,0,±1,±1}(r) = ±Li₃,₁-type, H_{0,±1,0,±1}(r) = ∓Li₂,₂-type (numerically verified). The explicit solve F = (Li₄ combination) + (Li₃,₁/Li₂,₂ combination) + (shuffle products), its 50‑digit verification at 3 angles, and the weight ≤ 3 residual identification were **not completed** in the available time; I do not claim any unprinted formula.

## 6. Task 5 — controls (partial)

Done (exact): Li₄ → 0; Li₃,₁(x,1) → 0; Li₂,₂(x,1) → 0; shuffles e₁шb³, (01)ш(10), Li₂(b)² → 0; positive control H_{0,1,0,−1}(−x) → nonzero; paper's full product symbol → 0. The scripted 100‑random‑control suite and the aggregate `verify_cor4.py` (< 5 min) were not completed.

## 7. Evidence table

| Check | Label | Result |
|---|---|---|
| step10: symbol conventions vs Li₄/shuffle/Li₃,₁/Li₂,₂ (exact) | PROVED (symbol) | all pass |
| step10: HPL series vs polylog/MPL, x=0.31,0.77 | NUMERICAL | max residual 1.4·10⁻⁶⁰ |
| S(F) exact, 56 terms (F_symbol.json) | PROVED (symbol) | printed in full |
| δ₂,₂(F) ≠ 0 as function of c; (3,1),(1,3) ≠ 0 | PROVED (symbol) | 2 components: −c⁻¹+1, −c⁻¹−1 |
| δ₂,₂ at c=1/√2 = 0 (required by thm:W4R) | PROVED (symbol) | 0 ✓ |
| δ₂,₂ at c=1/2, 1/3 (automatic, Λ²B₂=0) | PROVED (symbol) | 0 |
| δ₂,₂ at √3/2, 2/5, 3/7, 4/7, 5/9, 5/11, cos(3π/8), cos(π/12) | PROVED (symbol) | all ≠ 0 (values printed) |
| quartic letter relations 2=d²(2−d)(2+d), (2+d)(2−d)=d⁻² | PROVED (sympy exact) | == 0 |
| quartic basis independence | NUMERICAL | det ≈ −1.97; rank 3 |
| F(1/√2) vs paper W₄: D = explicit products | NUMERICAL | residual 3.61·10⁻⁵⁹ |
| M₃(G1) at 1/√2 vs paper M₃ | NUMERICAL | 2.61·10⁻⁵² |
| Li₄‑span membership at c=1/2 etc.; explicit depth‑2 identity + 50‑digit checks; 100 random controls; verify_cor4.py | — | not completed |

## 8. Files (in new_results/COR4_M3alpha_classical/)

`cor4_sym.py` (exact symbol algebra: HPL/Li symbols, canon, shuffles, δ₂,₂→Λ²B₂, δ₃,₁/δ₁,₃, ℚ(√2) S‑unit helper), `hpl_series.py` (fast high‑precision HPL/MPL evaluation), `step10_symbol_checks.py` (Task 1, PASS), `step11_F_symbol.py` (Task 2), `step12_specializations.py` (Task 3), `step13_library.py` (S‑unit args), `F_symbol.json`, `delta22_specializations.json`, `sunit_args.json`, plus auto‑logged `step_NN.py` transcripts of all computations.

**Bottom line:** M₃(α) is not classical as a function of the corner angle: the (2,2) Goncharov component of its weight‑4 part is the nonzero tensor (1−c⁻¹)([c]∧[1−c])∧([1−c]∧[1+c]) + (−1−c⁻¹)([c]∧[1+c])∧([1−c]∧[1+c]). The square corner α = π/2 is special because its letters collapse to {2, 1+√2}, forcing Λ²B₂ = 0 — the obstruction vanishes identically there, in agreement with the paper's proved classicality; at every other tested angle with ≥ 3 independent letters (π/3, 3π/4, π/6, and generic rational c) the obstruction persists.