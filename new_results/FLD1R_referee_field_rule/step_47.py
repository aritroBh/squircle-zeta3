result = r"""# FLD1R — adversarial referee report on FLD1's proof of the Field rule (conj:field) for every k

Referee workspace: `research/08_theory/FLD1R_referee_field_rule/`. All FLD1 sources were copied read-only and
re-run in this sandbox; all k=7 computations are the referee's own (independent re-implementations, cross-checked
against FLD1's code and against the raw corner integral). Every number below was printed by a script in this
directory (auto-archived `step_NN.py`; see the index at the bottom). `verify_fld1r.py` re-runs the essential
checks in **11 s (< 5 min), 5/5 PASS**.

## Verdict (one line)

**The Field rule as an existence theorem (HPL form with letters {0,±1} and denominators dividing
c^a(1−c)^b(1+c)^d) is PROVED for all k — after two minor repairs; the weight bound wt(z)+|w| ≤ k+1 is PROVED
for all k after one further (supplied) repair; the refinement "denominators are powers of c" and the sharp
pole set are PROVED (computationally) exactly for k ≤ 7, and are OPEN for k ≥ 8** (FLD1's own range was k ≤ 6;
this referee extended the exact verification to k = 7).

## Line-by-line lemma review (general k)

### Lemma A — exact nested-sum form of μ_{k,n}: **CORRECT (structure, A1) / MINOR FIX (A2 pole bound) / GAP (A2 properness)**
- (i) Beta-integral formula: correct for all k (∫ y^{n+j−1}(1−y)^{n−1} a^i b^ℓ = ∂_p^i∂_q^ℓ B|_{(n+j,n)};
  B(n+j,n)/B(n,n) = (n)_j/(2n)_j). Independently confirmed at k = 7: referee's own per-n series implementation
  of the Beta integrand == FLD1's Stirling implementation (exact, k ≤ 7, n ≤ 30, all 15 ζ-monomials), and the
  flat sum Σμ_{k,n}c^{2n} at c = 1/2 reproduces the **direct corner integral** `corner_series` M_k(2π/3) to
  ~22 digits (truncation-limited) for k = 5, 6, **7**.
- (ii) u-expansion: correct; referee notes the free sharpening **i+ℓ = k exactly** (checked k = 2, 3, 5):
  A^i ℓ^j has a,b-degree i+j = k, so D_{i,ℓ} only occurs at total Bell weight k.
- (iii) cumulants/polygammas: correct (exponential/Bell formula; ψ^{(m)}(N) = (−1)^{m+1}m!(ζ(m+1)−H^{(m+1)}_{N−1});
  pure cumulants ζ-free, mixed carry a single ζ(r+m), r+m ≥ 2 — hence only ζ(2),…,ζ(k), never ζ(1)).
- (iv)–(v) reduction to Z_w(2n) with |w| ≤ k−wt(z) (A1): correct for all k (stuffle preserves weight;
  H^{(r)}_{n−1} = 2^{r−1}(Z_r+Z_{−r})(2n) − n^{−r} and H^{(r)}_{2n−1} = Z_r(2n) − (2n)^{−r} verified exactly).
- (A2) **pole set — MINOR FIX.** The stated bound {1/2,0}∪{−t/2: 1 ≤ t ≤ k−1} is not what the given assembly
  proves: the shift corrections (n+t)^{−r}, 1 ≤ t ≤ j−1 ≤ k−1, give poles at n = −t = −(2t)/2 with 2t up to
  **2k−2**, not k−1. The argument proves the coarser (sufficient) set {1/2,0}∪{−t/2: 1 ≤ t ≤ 2k−2}; the sharp
  bound requires cancellations between (j,i,ℓ)-terms, not argued in the text. The sharp bound is however
  **verified exactly for k ≤ 7** (k = 7 poles: {−3,−5/2,−2,−3/2,−1,−1/2,0,1/2}). Either version suffices for
  every downstream claim: all poles lie in ½ℤ, hence are integral in m = 2n (what Lemma C needs).
- (A2) **properness (deg num ≤ deg den −2) — GAP for general k** (PROVED computationally for k ≤ 7). The
  offered argument ("the analytic statement μ_{k,n} = O(n^{−2})") is insufficient as stated: Z_w(2n) grow only
  poly-logarithmically, so decay of the total sum does not imply per-word properness without a
  **linear-independence/separation lemma** for the sequences {Z_w(2n)} over ℚ(n) (true — e.g. via Schneider's
  difference-ring theory of alternating harmonic sums — but neither cited nor proved by FLD1). Exact check:
  worst deg = −2 for every ρ_w at k = 7 (551 word-terms). This gap affects ONLY the "powers of c denominators"
  refinement — NOT existence (Lemma C's polynomial-part clause covers it).

### Lemma B — the dictionary w ↔ (ε,k,v): **CORRECT (all lengths)**
Induction/recursion checked line by line (base h_(a)(N) = a^{N+1}/N in the text should read
h_(a)(N) = a^{N−1}/N — a typo only, cf = a, ε = a, k = 1, v = () is what is used and is correct;
h_{0,w} = h_w/N; h_{a,w}(N) = (1/N)Σ a^{N−1−i}h_w(i)). Injectivity/surjectivity argument is complete.
Independently re-verified against **direct Taylor coefficients** (hpl.coeffs_exact, which implements the HPL
integral recursion, not the formula): all **728 words of length ≤ 6**, N = 3..13; map injective, inverse
reconstructs w, |w| = k+|v|.

### Lemma C — the construction: **CORRECT for all k, with two MINOR FIXES**
- Termination: correct (each step cancels a longest word exactly and adds only strict suffixes; the
  suffix-closure of a finite word set is finite; the killed word can never be re-created).
- The bijection term c^e H_{bji(1,κ,w)} has leading sum cf·(2n−e)^{−κ}Z_w(2n) exactly (hpl_coeff_at + _shift,
  whose corrections are strict suffixes) ✓.
- MINOR FIX 1: the parenthetical "no denominator 2n+t used in the shifts vanishes there [n ≥ 1]" is false as
  stated (negative shifts produce 2n−t, t ≤ e, which vanishes at n = t/2 — e.g. n = 1 when t = 2). Repair:
  the identity is an identity of rational functions; the Z-sum expression E(n) is regular at every integer
  n ≥ 1 (its poles are at {1/2,0} ∪ negative half-integers) and [c^{2n}]F is a plain rational number, so the
  identity — valid at all but finitely many n — extends to every n ≥ 1. Conclusion unchanged.
- MINOR FIX 2: polynomial parts P(n) are handled in the text by (c d/dc)^j but are **not implemented** in
  step_25_construct.py (it raises on polynomial parts). The textual argument is sound ((c d/dc)^j applied to
  the whole realizing combination multiplies [c^{2n}] by (2n)^j; c d/dc preserves the class
  {r(c)H_w : den | c(1−c)(1+c)} and does not increase word length). Vacuous for k ≤ 7 (no polynomial parts occur).
- Symmetrization F ↦ (F(c)+F(−c))/2 with H_w(−c) = (−1)^{nz(w)}H_{ν(w)}(c): correct.
- **Finite low-order correction: CORRECT.** Coefficients of F_1 at nonpositive even powers come only from the
  finitely many terms with prefactor c^e, e ≤ 0; for j ≤ min(0,e), [c^j]c^eH_w is the Laurent tail of the
  rational product of the kernels — a Laurent polynomial; the sum is a finite even Laurent polynomial L(c)
  (deg bounded by e_min ≥ −(k−1) using the sharp pole set; ≥ −(2k−2) with the coarse one); odd negative
  coefficients vanish since F_1 is even. Subtracting L fixes c^0 via M_k(0) = 0. Verified at k = 7
  (corrections: rational part none; ζ(2)-part c^0: −1/2; ζ(3)-part c^0: +1/2 — matching FLD1's "∓1/2"
  pattern; single-ζ(6,7) parts need c^{−4}, c^{−2}).
- **Convergence and identity of analytic functions: CORRECT.** Each H_w (last letter ≠ 0) is analytic on
  |c| < 1 with H_w(0) = 0; after Step 3 the total F has no negative powers, so F is analytic on |c| < 1; the
  μ-side converges there by thm:flat; coefficientwise equality (all N, via the repaired Step-1 argument) gives
  F ≡ M_k on |c| < 1 ⊃ (0,1) ∋ c = cos(α/2) for 0 < α < π.

### Lemma D — weight bound wt(z)+|w_HPL| ≤ k+1: **PROVED AFTER REPAIRS (all k)** (upgrade of FLD1's "k ≤ 6")
The budget argument as written **forgets the simple poles of the ratio (n)_j/(2n)_j**. Repair (supplied by the
referee, checked): (n)_j/(2n)_j = Π_{t=0}^{j−1}(n+t)/(2n+t) has NO pole at 0 (the t = 0 factor is n/(2n) = 1/2)
and simple poles only at −t/2, 1 ≤ t ≤ j−1 — disjoint from the prefactor's poles {0, 1/2}. Hence at every
pole location n₀, ord_{n₀}(ρ_w) ≤ (Σ shift-orders) + 1 ≤ (k − wt(z) − |w|) + 1, so
|w| + ν(ρ_w) ≤ k+1−wt(z); the HPL word w0 = bji(1,κ,w) has |w0| = κ+|w| ≤ ν+|w| ≤ k+1−wt(z) (Lemma B2);
polynomial parts, if any, are realized by c d/dc, which does not increase word length. Therefore
wt(z)+|w_HPL| ≤ k+1 for ALL k. Verified per word AND per pole location at k = 7: 0 violations, saturated
(max = bound for every ζ-monomial).

### Lemma E — letters {0,±1}, −1 necessary: **CORRECT (all k)**
bji outputs {0,±1}-words, ν permutes them, c d/dc adds none. Necessity of −1 via the even/odd split
H^{(r)}_{n−1} = 2^{r−1}(Z_r+Z_{−r}) − n^{−r} is correct; at k = 7 the rational part uses −1 in 4513/4678
(96.5%) of terms.

### Lemma F — S-units ⟨2,c,1−c,1+c⟩: **CORRECT**
H_w(c) are iterated integrals of dlog t, dlog(1−t), dlog(1+t) at t = c; as MPLs the arguments are ±c, ±1
(trivially S-units), and the symbol letters are exactly {2, c, 1−c, 1+c} (checked with hpl.symbol). This is
conj:field with the explicit S = {primes of ℚ(cos(α/2)) dividing 2, c, 1−c, 1+c}.

## Task 2 — independent test at k = 7 (beyond FLD1's k ≤ 6)

| Check | Result |
|---|---|
| μ_{7,n} exact, n ≤ 30, referee's per-n sympy-series implementation vs FLD1's Stirling implementation | **identical** (all 15 ζ-monomials; `step_26`, `step_28`) |
| Σ_{n≤30} μ_{k,n} c^{2n} at c = 1/2 vs **direct corner integral** (corner_series, k = 5,6,7) | residuals 5.2e-23 / 4.6e-23 / **3.9e-23** (~22 digits, truncation-limited) (`step_30`) |
| symbolic-n μ_formula(7) (429 monomials) vs exact table | **identical**, n ≤ 30 (`step_31`) |
| Z-sum form k = 7 (551 word-terms): ρ_w proper | **all deg = −2 exactly** (`step_34`) |
| poles of ρ_w at k = 7 | {−3,−5/2,−2,−3/2,−1,−1/2,0,1/2} = sharp set {1/2,0}∪{−t/2: t≤6} (`step_34`) |
| (A1) |w| ≤ 7−wt(z) | 0 violations (`step_34`) |
| budget |w|+ord ≤ 8−wt(z) per word & pole location | 0 violations, **saturated** (`step_34`) |
| construction at k = 7 (all 15 monomials): all-n even identity == μ AND odd == 0 (ℚ(n) identities) | **True** for all 15 (`step_35`, `step_36`) |
| exact [c^N]F == μ_{7,n}, N = −6..60 (n ≤ 30), odd/nonpositive == 0 | **True**, all 15 monomials (`step_38`) |
| final HPL form (6553 terms): letters, denominators, weight bound | letters {−1,0,1}; min prefactor e = −6 (powers of c only); max wt(z)+|w| = 8, 0 violations (`step_37`) |
| k = 7 HPL form vs direct corner integral M_7(2π/3) | residual **3.5e-46** (~45 digits) (`step_39`) |
| FLD1's own verify_fld1.py rerun in this sandbox | **5/5 PASS in 74 s** (`step_40`) |
| FLD1 k = 3 Z-sum form vs W1 template mu_zsum.pkl | identical, all 3 monomials, all words (`step_41/42`) |

## Task 3 — is the THEOREM (existence for all k) proved, with weight bound/properness only k ≤ 6?

**Yes for existence** (now even with the weight bound for all k, after the Lemma-D repair): Lemmas A (structure),
B, C, E, F are correct for general k, and Lemma C's polynomial-part clause makes the construction independent
of the unproven properness — at the price of (1±c) denominators should a polynomial part ever occur. **No for
the exact statement FLD1 printed** ("denominators divide a power of c" for all k): that clause is properness
(deg ≤ −1 suffices, deg ≤ −2 verified), proved computationally only for k ≤ 7. Note FLD1's labels were overly
conservative for Lemma D (the budget argument IS general once the ratio poles are accounted) and overly
optimistic for the (A2) pole bound (the text's assembly gives t ≤ 2k−2, not t ≤ k−1).

## Task 4 — verdict table and proposed wording

| Claim | Verdict |
|---|---|
| Existence of the HPL form, letters {0,±1}, denominators | c^a(1−c)^b(1+c)^d, **all k** | **PROVED** (after the two Lemma-C minor fixes) |
| Weight bound wt(z)+|w| ≤ k+1, **all k** | **PROVED AFTER REPAIRS** (ratio-pole accounting; verified k ≤ 7, saturated) |
| Denominators are powers of c (properness, deg ≤ −2) | **PROVED for k ≤ 7** (exact computation); **NOT PROVED for k ≥ 8** (needs a Z-sum separation lemma) |
| Sharp pole set {1/2,0}∪{−t/2: t≤k−1} | **PROVED for k ≤ 7**; for general k only the coarser t ≤ 2k−2 follows from the text (sufficient) |
| conj:field with S = {2, c, 1−c, 1+c}, all k | **PROVED** (corollary of the first row) |

**Proposed Theorem wording for the paper.**

> **Theorem (Field rule).** For every integer k ≥ 0 and every 0 < α < π, with c = cos(α/2),
> M_k(α) = Σ_j z_j r_j(c) H_{w_j}(c), a finite sum in which (i) the z_j are monomials in ζ(2),…,ζ(k);
> (ii) the r_j ∈ ℚ(c) have denominators dividing a product c^a(1−c)^b(1+c)^d; (iii) the H_{w_j} are harmonic
> polylogarithms with letters in {0,1,−1}; (iv) weight(z_j) + |w_j| ≤ k+1. Moreover, for every k ≤ 7 the
> denominators divide a power of c alone.
>
> **Corollary (conj:field).** Every M_k(α) is a combination of multiple polylogarithms at S-units of
> ℚ(cos(α/2)), with S = {primes dividing 2, c, 1−c, 1+c}.
>
> Remark. The power-of-c denominators and the pole set {1/2,0}∪{−t/2: 1≤t≤k−1} hold for all k ≤ 7 by exact
> computation and are conjectured for all k; they follow from the properness of the coefficients ρ_w(n) in
> Lemma A, which we verify for k ≤ 7.

## Files produced by this referee (work dir)
- `verify_fld1r.py` — 5/5 PASS in 11 s (< 5 min).
- Archived steps: `step_23.py` (setup + zsum self-tests), `step_24.py` (Lemma B independent, 728 words),
  `step_26.py` (independent μ table k ≤ 7, n ≤ 30 → `mu7_indep_N30.json`; note `step_25.py` was a timed-out
  first attempt with symbolic n), `step_28.py` (Stirling cross-check → `mu_exact_K7_N30.json`),
  `step_30.py` (corner_series k=5,6,7), `step_31.py` (symbolic μ_7(n) → `mu_symbolic_K7.pkl`),
  `step_32.py` (Z-sum form → `mu_zsum_K7.pkl`), `step_34.py` (properness/poles/budget),
  `step_35.py` + `step_36.py` (construction k = 7 → `hpl_constructed_K7.pkl`),
  `step_37.py`/`step_38.py` (weight bound + exact coefficients), `step_39.py` (k=7 end-to-end, 45 digits),
  `step_40.py` (FLD1 verify_fld1.py rerun), `step_41.py`/`step_42.py` (W1 template + Lemma A(ii) support),
  `step_43.py`/`step_44.py` (verify_fld1r write/run), `step_45.py` (summary stats).
- Artifacts: `mu7_indep_N30.json`, `mu_exact_K7_N30.json`, `mu_symbolic_K7.pkl`, `mu_zsum_K7.pkl`,
  `hpl_constructed_K7.pkl`, `hpl_constructed_K7_nonrat.pkl`, `progress.log`.
"""
open('result.md','w').write(result)
open('progress.log','a').write('step09: result.md written (full referee report)\n')
print('result.md written:', len(result), 'chars')
