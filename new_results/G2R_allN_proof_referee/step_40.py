CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
md = r'''# G2R referee report — the all-n proof of Theorem B4 (M₄ at every corner angle)

Work dir: `research/08_theory/G2R_allN_proof_referee`. Refereed object: `G2_M4alpha/step_11.py`
(with `step_10.py`, `g2_fit.py`, the term lists `hpl_fit_k4_(0,0,0).json`, `hpl_fit_k4_(1,0,0).json`,
G1's `hpl_fit_k3_(1,0,0).json`, `b4_pi4_part.json`, `mu4_formula.pkl`), following the Z-sum method of
`W1_corner_family/B_taylor_pi4/b3_allN_proof.py` (`zsum.py`). Every number below was printed by a
script in this folder (see §7). **Overall verdict: PROVED** (no repairs needed).

## 1. What step_11.py proves (method, in my own words)

For each zeta sector P ∈ {1, ζ(2), ζ(3), ζ(4), ζ(2)²} separately, let B4|_P be the claimed HPL
combination (R₄; Ψ₂; −Ψ with Ψ = B3's ζ(2)-part; F₁; F₂ — all given as finite sums of
q·c^e·H_w(c), plus the even-sector reading He[w] := c^e(H_w + (−1)^{e+nz(w)}H_{−w}), plus F₁,F₂'s
Laurent monomials) and let μ₄,ₙ|_P be the matching part of step_05's symbolic μ₄,ₙ. step_11 proves,
for every integer n ≥ 1:

  (even)  [c^{2n}]  B4|_P  =  μ₄,ₙ|_P ,
  (odd)   [c^{2n+1}] B4|_P =  0 ,

and it checks the remaining coefficients c^{−3}…c³ directly in exact arithmetic (targets: 0 except
[c²] = μ₄,₁|_P). Coverage is complete: even N ≥ 2 and odd N ≥ 3 by the symbolic identities,
N ∈ {−3,…,1} by the direct check (N = 2, 3 redundantly both ways); no term has valuation < −3
(min e = −3 with H_w = O(c), F₁/F₂ monomials stop at c^{−2}), so every coefficient is accounted for.
With Theorem A (M₄(α) = Σ_{n≥1} μ₄,ₙ c^{2n}, |c| < 1) and polynomial boundedness of μ₄,ₙ, both sides
are analytic on |c| < 1 (HPLs/artanh have their first branch points at c = ±1; the Laurent parts
cancel by the coefficient check), so equal coefficients ⇒ B4 = M₄ as functions. The per-part split is
*stronger* than needed: no ℚ-linear-independence assumption on {1, ζ2, ζ3, ζ4, ζ2²} is used, since
each part-difference is separately the zero function.

**Which Z-sums, which upper limits.** The μ side is a ℚ(n)-polynomial in H⁽ʳ⁾_{n−1} and
H⁽ʳ⁾_{2n−1} (r ≤ 4) — upper limits n−1 and 2n−1. The HPL coefficient formula is
h_w(N) := [c^N]H_w = Σ coef·ε^N N^{−k} Z_v(N−1), i.e. Z-sums at N−1 = 2n+δ−1 with
δ = parity − e ∈ {−1,…,4}. Everything is canonicalised to alternating Z-sums at the single upper
limit 2n: Z_w(2n), w a word over ℤ\{0} (observed indices {−4,…,−1, 2,…,4}, depth ≤ 4; canonical
forms have 29/5/2/1/1 words for the five parts; the odd forms cancel to the empty combination).

**Which identities the canonicalisation uses** (all exact, finite rearrangements, hence valid for
every n where defined):
  1. quasi-shuffle (stuffle): Z_u(M)Z_v(M) = Σ_w Z_w(M) with merge a∘b = sgn(a)sgn(b)(|a|+|b|) —
     used to expand products of harmonic sums in μ₄,ₙ;
  2. even/odd index splitting: H⁽ʳ⁾_{n−1} = 2^{r−1}(Z_r(2n) + Z_{−r}(2n)) − n^{−r} (the even terms
     of H⁽ʳ⁾_{2n} give H⁽ʳ⁾_n), and H⁽ʳ⁾_{2n−1} = Z_r(2n) − (2n)^{−r};
  3. index shifts: Z_{a,v}(2n+j) = Z_{a,v}(2n) + Σ_t sgn(a)^t (2n+t)^{−|a|} Z_v(2n+t−1), both signs
     of j, applied recursively (j ∈ {−2,…,3} occurs);
  4. the HPL Taylor recursion from the integral definition: h_{±1}(N) = (±1)^{N∓1}/N,
     h_{0,w}(N) = h_w(N)/N, h_{a,w}(N) = (1/N)Σ_{i<N} a^{N−1−i} h_w(i) = coef·a^{N±1}/N·Z_{(aεk),v}(N−1).
  Partial fractions in the summation index are *not* used inside zsum; they appear only upstream in
  step_07's derivation of F₁, F₂ (apart decomposition of (n+3)/((2n+3)(4n²−1)) etc. into
  1/(2n+m) generators with artanh closed forms) — re-verified by me from the pickle (identities hold;
  apart forms match).

**Validity for every n ≥ 1, boundaries included.** The shifted sums keep upper limits
2n+j ≥ 0 for n ≥ 1, j ≥ −2 (Z(0) = 0 for nonempty words, exact); all denominators ever introduced
are factors of {n, n+1, 2n−1, 2n+1, 2n+3} (my audit over all ten canonical forms found exactly
these), none vanishing at an integer n ≥ 1. I additionally verified *every* symbolic-identity
instance used, exactly, at the boundary: hpl_coeff_at(w, δ) vs the direct Taylor coefficient for all
96 words of B4 (incl. negations), δ ∈ {−1,…,4}, n = 1,2,3 — 1728 checks, all True; shift(v, j) vs
Z_v(2n+j) for all 29 Z-words used, j ∈ {−2,…,3}, n = 1,2 — 348 checks, all True; harmonic
conversions r ≤ 4 at n = 1,2,3 — True.

**"Identical canonical forms ⇒ equality"**: yes, and only this direction is needed. The check is
`diff == {}` — the difference of the two sides simplifies to the *literally zero* ℚ(n)-combination
(0 coefficient of every word). Since each rewrite above is an exact identity for every n ≥ 1, a
zero combination evaluates to zero for every n ≥ 1. Linear relations among Z-sums (which exist)
only affect uniqueness of the form, never the sufficiency. Coefficient arithmetic is exact in the
sympy rational-function field ℚ(n).

**HPL → Z-sum expansion and the c^e prefactors (e ∈ {−3,…,1})**: used correctly.
[c^N]c^eH_w = h_w(N−e) with N = 2n (even) resp. 2n+1 (odd); negative e shifts the index to
N−e > N (δ down to −1 for even, up to 4 for odd) — this is exactly what hpl_coeff_at(w, parity−e)
computes, with ε^{2n+δ} = ε^δ (sign bookkeeping via `eps**(delta%2)`), and the even-sector partner
H_{−w} enters with sign (−1)^{e+nz(w)}, consistent with g2_fit's `even_series` and with the
planner's independent reading (H_w(−c) = (−1)^{nz(w)}H_{−w}(c)). Cross-validated: my own exact HPL
evaluator (written from the integral definition) agrees with `hpl.coeffs_exact` for all 96 words,
m < 84, and both agree with `hpl_coeff_at` evaluations at n = 1,2,3 (above).

**Low orders c^{−3}…c³**: checked by `low_order_check` per part in exact arithmetic (kinds
'1'/'c'/'1c'/hpl/even all handled; hpl terms contribute only when N−e ≥ 1); targets 0 except
[c²] = μ₄,₁|_P from the table. I re-verified N = −3…1 against the *definition* independently
(spot check 4c).

## 2. The run (task 2)

Pristine copy in `run_base/`, run with the project venv: **rc = 0, wall 3.7 s** (script self-report
2–3 s). Output (`run_base.out`): for each part, `even identity … True`, `odd identity … True`,
`low orders … True`; final line `B4 PROVED for all n (all five parts): True`. Canonical-form sizes
29/5/2/1/1 as claimed. The regenerated `mu4_zsum.pkl` equals the stored artifact exactly.

**Where equality is asserted** (`step_11.py` __main__): `ev = hpl_side(terms, 0)` (built from the
JSON term lists = the B4 candidate), `muE` from `Z = mu4_zsum()` (built from `mu4_formula.pkl` =
step_05's independent symbolic-n evaluation of the μ-integral); then `diff = zsum.lin((1, ev),
(-1, muE))` must satisfy `len(diff) == 0`; `od = hpl_side(terms, 1)` must satisfy `len(od) == 0`;
`lo = low_order_check(...)`; `allok = allok and not diff and not od and lo`.

**Not vacuous.** (i) The two sides come from independent inputs (JSON fits vs the μ formula) and
both are nonempty (29/5/2/1/1 words). (ii) The odd identity of part '1' alone accumulates **446
(word, coefficient) insertions that cancel to 0 surviving words**; the even side at n = 5 is an
exact rational with L1-magnitude 3.53·10⁻¹ cancelling to the value 3.46·10⁻², agreeing with the μ
side to the exact fraction −7596252093627541/219583673971200000 (= flat_exact's μ_{4,5} rational
part — three-way agreement). (iii) The mutation tests below all fail the proof — the comparison
reacts to any change of input.

## 3. Mutation tests (task 3)

Each mutant = full copy of the pipeline with one input perturbed; ran `step_11.py` (rc = 0 each,
~3.5 s). **All four FAILED, as a correct proof must.** In each case only the affected part breaks
(the other four still pass — per-part sensitivity), and the failure is caught both symbolically
(nonempty residual) and by the low-order check:

| mutant | perturbation | even residual (word: ℚ(n) coef) | low-order mismatch | verdict |
|---|---|---|---|---|
| (a) R₄ weight-5 coefficient | `c^0 He[0,0,0,0,1]`: 1 → 2 | Z_() : 1/(16n⁵) | c²: −1007/2160 vs −571/1080 | FAILED ✓ |
| (b) Ψ₂ coefficient | `c^0 He[1]`: 5/4 → 19/12 | Z_() : 1/(3n) | c²: 247/432 vs 103/432 | FAILED ✓ |
| (c) F₁ coefficient | artanh prefactor 7c/16 → 5c/16 | Z_() : −1/(16n−8) | c²: 17/120 vs 4/15 | FAILED ✓ |
| (d) prefactor exponent | `c^-2 He[1,0,1]` → `c^-3 He[1,0,1]` (coef 2) | Z_() : (4n²+2n−1)/(8n⁵+36n⁴+62n³+51n²+20n+3); Z_(2) : −2/(2n²+5n+3) | c⁰: −1/3 vs 0; c²: −811/1080 vs −571/1080 | FAILED ✓ |

The odd identity stays True in all four mutants, as it should: the He basis functions are even by
construction, and the mutated F₁ remains even (odd×artanh). No mutant passes ⇒ the proof is not
vacuously true, and it detects coefficient-level and exponent-level changes alike.

## 4. Independent spot check (task 4)

`spot_check.py`, fully independent code (own parser, own exact HPL Taylor evaluator from the
integral definition; μ₄,ₙ from the **definition** via G1's unmodified `flat_exact.py 40 4`, run by
me in `flat_run/`, 92 s):

- (4a) [c^{2n}](B4 part) == μ₄,ₙ part **exactly in ℚ, per zeta monomial, all n = 1…40, all five
  parts: True**;
- (4b) odd coefficients c¹…c⁷⁹ vanish, all parts: True;
- (4c) low orders c^{−3}…c¹ vanish, all parts: True;
- (4d) the proof's μ input (`mu4_formula.pkl`, step_05's symbolic formula) evaluated at the
  harmonic numbers == the flat_exact definition for all n ≤ 40, all parts: True (so the object the
  proof compares against *is* μ₄,ₙ).

Corroboration reproduced in this folder: `step_10.py` exact check **EXACT through c^1280, all five
parts: True** (30 s); `zsum.py` self-tests True; planner's independent
`planner_check_b4.py`: |B4 − M4_800| = **8.6732e-131** at c = 1/√2 (dps 130), 4.8651e-93 at 2π/3,
2.8715e-92 at π/3; step_07's [ζ(4)]μ = (n+3)/((2n+3)(4n²−1)) and [ζ(2)²]μ = −n/((2n+3)(4n²−1))
re-derived from the pickle: both True.

## 5. Verdict table

| item | verdict | evidence |
|---|---|---|
| Method: upper limits n−1, 2n−1 → canonical Z(2n) | CORRECT | §1; referee_checks (2a)(2b)(2c): 1728+348+harmonic checks at n=1,2,3, all True |
| Method: stuffle / even-odd splitting / shifts / HPL coefficient recursion | CORRECT | exact identities; boundary-verified; zsum self-tests True |
| Validity every n ≥ 1 (denominators, ranges, small n) | CORRECT | denominator audit: only {n, n+1, 2n−1, 2n+1, 2n+3}, no integer zeros ≥ 1; n=1,2,3 boundary checks True |
| "Identical canonical forms ⇒ equality" | CORRECT | diff is the zero combination; only exact rewrites used; reverse direction not needed |
| c^e prefactors e ∈ {−3,…,1} (index shift, ε^δ signs, He sign) | CORRECT | §1; own evaluator == coeffs_exact == hpl_coeff_at, 96 words |
| Low orders c^{−3}…c³ | CORRECT | low_order_check True in the run; independently re-checked (4c) |
| The run (rc, time, verdict) | REPRODUCED | rc 0, 3.7 s, `PROVED … True`, run_base.out |
| Non-vacuity | CONFIRMED | 446 insertions → 0; three-way value agreement at n=5; mutants all detected |
| Mutant (a) R₄ w5 coefficient | FAILED as required | residual 1/(16n⁵), low-order mismatch |
| Mutant (b) Ψ₂ coefficient | FAILED as required | residual 1/(3n), low-order mismatch |
| Mutant (c) F₁ coefficient | FAILED as required | residual −1/(16n−8), low-order mismatch |
| Mutant (d) prefactor exponent | FAILED as required | residuals on Z_() and Z_(2), low-order mismatches at c⁰, c² |
| Spot check n = 1…40, exact, per part | PASS | (4a)–(4d) all True |
| G2 result.md §2 weight census "3/8/28/70" | COSMETIC GAP (documentation only) | actual JSON census is 1/8/27/73 (total 109 either way); does not touch the proof |

## 6. Overall verdict: **PROVED**

The all-n proof of Theorem B4 is correct, complete on every coefficient, non-vacuous, and robust
under mutation. Scope notes (not gaps of step_11): (i) the proof takes as inputs Theorem A
(M₄ = Σ μ₄,ₙ c^{2n} on |c| < 1) and step_05's exact symbolic-n evaluation of μ₄,ₙ — the latter is
itself cross-validated against G1's independent `flat_exact.py` for all n ≤ 640 (step_05) and
n ≤ 40 by me (4d); (ii) G2's result.md §2 misstates the weight census of R₄ (documentation only).

## 7. Files in this folder

- `verify_g2r.py` (< 5 min; actual 7 s): pristine pipeline copy → run `step_11.py` (must say
  `True`), then mutant (a) (must say `False`); exit 0 iff both behave. `verify_g2r.out`: PASS.
- `referee_checks.py` / `.out`: non-vacuity (446 insertions), boundary checks (1728 + 348),
  denominator audit, own-evaluator cross-check (96 words × m < 84).
- `spot_check.py` / `.out`: task-4 exact spot check n = 1…40 + (4d); `flat_run/` holds the
  `flat_exact.py 40 4` run (the definition of μ₄,ₙ).
- `run_base/` (pristine pipeline copy) with `run_base.out` (step_11) and `run_base_step10.out`
  (EXACT-1280); `mut_a/`…`mut_d/` with `mut_a.out`…`mut_d.out` (the four mutants).
- `progress.log`: one line per finished step.
'''
open(CWD + "/result.md", "w").write(md)
with open(CWD + "/progress.log", "a") as fh:
    fh.write("step_result: result.md written; verdict PROVED\n")
print(open(CWD + "/result.md").read()[:600])
print("...")
print(open(CWD + "/progress.log").read())
