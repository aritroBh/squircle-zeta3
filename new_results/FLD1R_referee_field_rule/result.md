# FLD1R referee verdict on FLD1's proof of the Field rule (conj:field) for every k

## Verdict table

| Item | Verdict |
|---|---|
| **Lemma A** (Beta-integral structure; (A1) |w| ≤ k−wt(z)) | **CORRECT** for all k (verified exactly k ≤ 7; free sharpening: i+ℓ = k *exactly*) |
| Lemma A (A2), pole set {1/2,0}∪{−t/2: t≤k−1} | **MINOR FIX**: the text's assembly only proves the coarser set {1/2,0}∪{−t/2: t ≤ 2k−2} (the (n+t)-shifts give −t = −(2t)/2, 2t ≤ 2k−2); the sharp bound needs un-argued cancellations. Sharp set verified exactly for k ≤ 7. Either version suffices downstream (all poles in ½ℤ ⟹ integral in m=2n). |
| Lemma A (A2), properness deg ≤ −2 | **GAP for general k**: "μ_{k,n}=O(n^{−2}) analytically" is insufficient — per-word properness needs a linear-independence/separation lemma for {Z_w(2n)} over ℚ(n) (available in principle via Schneider's difference-ring theory, but not cited or proved). **PROVED computationally for k ≤ 7** (worst deg = −2, all 553 word-terms at k=7). Affects only the "powers of c denominators" refinement, not existence. |
| **Lemma B** (bijection w↔(ε,k,v), all lengths) | **CORRECT** — independently re-verified against direct Taylor coefficients (not the same recursion): all 728 words of length ≤ 6, N=3..13; injective, inverse reconstructs w, |w|=k+|v|. |
| **Lemma C** (construction terminates, exact, all k) | **CORRECT with two MINOR FIXES**: (1) the parenthetical "no denominator 2n+t vanishes at n≥1" is false (2n−2 vanishes at n=1); repair: the identity is an identity of rational functions, regular at every integer n≥1, hence extends there. (2) The polynomial-part clause (c·d/dc) is sound but unimplemented; vacuous for k≤7 (no polynomial parts). |
| Finite low-order correction | **CORRECT** (finite even Laurent polynomial from the ≤ 0-prefactor terms; odd negative parts vanish since F₁ is even; verified at k=7: rational part needs none; ζ(2): −1/2, ζ(3): +1/2 constant corrections). |
| Convergence for |c|<1; identity of analytic functions | **CORRECT** (F analytic on |c|<1 after correction; coefficientwise identity at all N ⇒ F≡M_k on |c|<1 ⊃ (0,1)). |
| **Lemma D** (weight bound wt(z)+|w|≤k+1) | **PROVED AFTER REPAIRS for all k** (upgrade of FLD1's "k≤6"): the budget argument forgets the simple poles of (n)_j/(2n)_j. Repair: the ratio has no pole at 0 (n/(2n)=1/2) and poles only at −t/2 (t≥1), disjoint from the prefactor's {0,1/2}; hence ord ≤ (k−wt(z)−|w|)+1 at every location, so |w|+ν ≤ k+1−wt(z), and |w_HPL|=κ+|w|≤ν+|w| (Lemma B2); c·d/dc never raises word length. Verified per word AND per pole location at k=7: 0 violations, saturated. |
| **Lemma E** (letters {0,±1}, −1 necessary) | **CORRECT** (all k; at k=7 the rational part uses −1 in 4513/4678 = 96.5% of terms). |
| **Lemma F** (S-units ⟨2,c,1−c,1+c⟩) | **CORRECT** (MPL arguments ±c,±1; symbol letters exactly {2,c,1−c,1+c}). |

## Answers to the referee questions

- **Q3 — is the THEOREM (existence for all k) proved even though weight bound/properness were only k≤6?** **Yes for existence**: Lemmas A, B, C, E, F are general-k correct and Lemma C's polynomial-part clause makes existence independent of the unproven properness (costing at most (1±c) denominators). Moreover, after the supplied Lemma-D repair, the **weight bound wt(z)+|w| ≤ k+1 is itself proved for all k**. **No for FLD1's printed statement** "denominators divide a power of c for all k": that clause is properness, proved computationally only for k ≤ 7 (FLD1 had k ≤ 6; I extended the exact verification to k = 7). Note FLD1's labels were overly conservative for Lemma D (the budget argument *is* general once the ratio poles are accounted) and overly optimistic for the (A2) pole bound.

- **Verdicts**: all-k HPL existence (denominators | c^a(1−c)^b(1+c)^d), letters, S-units, weight bound: **PROVED** (after the stated minor repairs). Power-of-c denominators and sharp pole set: **PROVED for k ≤ 7, NOT PROVED for k ≥ 8**.

## Proposed Theorem wording for the paper

> **Theorem (Field rule).** For every integer k ≥ 0 and every 0 < α < π, with c = cos(α/2), M_k(α) = Σ_j z_j r_j(c) H_{w_j}(c), a finite sum in which (i) the z_j are monomials in ζ(2),…,ζ(k); (ii) the r_j ∈ ℚ(c) have denominators dividing a product c^a(1−c)^b(1+c)^d; (iii) the H_{w_j} are harmonic polylogarithms with letters in {0,1,−1}; (iv) weight(z_j)+|w_j| ≤ k+1. Moreover, for every k ≤ 7 the denominators divide a power of c alone.
>
> **Corollary (conj:field).** Every M_k(α) is a combination of multiple polylogarithms at S-units of ℚ(cos(α/2)) with S = {primes dividing 2, c, 1−c, 1+c}.
>
> *Remark.* The power-of-c denominators and the sharp pole set {1/2,0}∪{−t/2: 1≤t≤k−1} are exact computations for k ≤ 7 and conjectured for all k; they would follow in general from the properness of the coefficients ρ_w(n) in Lemma A.

## Evidence table (all numbers printed by code; residuals and digit counts)

| Check | Result | Script |
|---|---|---|
| zsum self-tests (shift/stuffle/harmonic/HPL-coeff vs direct definitions) | True | step_23.py |
| Lemma B independent (728 words len ≤ 6, direct Taylor ground truth) | map/inverse/|w|=k+|v| all True | step_24.py |
| μ_{7,n} exact n ≤ 30: referee's per-n sympy-series impl. == FLD1's Stirling impl. | **identical**, all 15 ζ-monomials (μ_{7,1} printed, 15 keys) | step_26.py, step_28.py |
| Σ_{n≤30}μ_{k,n}c^{2n} at c=1/2 vs **direct corner integral**, k=5,6,7 | residuals 5.2e-23 / 4.6e-23 / **3.9e-23** (~22 digits, truncation-limited) | step_30.py |
| symbolic-n μ_formula(7) (429 monomials) == exact table, n ≤ 30 | True | step_31.py |
| k=7 Z-sum form (553 word-terms): proper; poles; (A1); budget | all deg = **−2**; poles {−3,−5/2,−2,−3/2,−1,−1/2,0,1/2}; 0 violations; budget saturated | step_32.py, step_34.py |
| k=7 construction, all 15 ζ-monomials: even==μ AND odd==0 (ℚ(n) identities) | True (rational part 4678 terms, correction {}) | step_35.py, step_36.py |
| exact [c^N]F == μ_{7,n}, N=−6..60; odd/nonpositive == 0 | True, all 15 monomials | step_38.py |
| k=7 HPL form (6553 terms): letters; denominators; weight bound | letters {−1,0,1}; min prefactor e=−6 (powers of c only); max wt(z)+|w|=8, 0 violations | step_37.py |
| k=7 HPL form vs direct corner integral M_7(2π/3), dps=50, N=500 | residual **3.5e-46 (~45 digits)** | step_39.py |
| FLD1's own verify_fld1.py rerun in referee sandbox | **5/5 PASS in 74 s** | step_40.py |
| FLD1 k=3 Z-sum form vs W1 template mu_zsum.pkl | identical (all 3 monomials, all words); Lemma A(ii) support: i+ℓ=k exactly (k=2,3,5) | step_41.py, step_42.py |
| **`verify_fld1r.py` (deliverable)** | **5/5 PASS in 10 s (< 5 min)** | verify_fld1r.py, step_43.py, step_44.py |

Artifacts: `mu7_indep_N30.json`, `mu_exact_K7_N30.json`, `mu_symbolic_K7.pkl`, `mu_zsum_K7.pkl`, `hpl_constructed_K7.pkl`, `result.md` (full report), `progress.log` (one line per step).