> **PLANNER CORRECTION (2026-10-22 19:08, from the CX4 Codex review, planner-verified):** Lemma B's displayed recursion is wrong for
> adjacent keys j = i−1: the τ2 feed must carry [i ≥ j+2], not [i > j], otherwise the diagonal feed is counted twice. Correct form:
> h_ij = δ_{k,d}ρ_ij + [(i+1)τ1_0 h_{i+1,j} + [i ≥ j+2](j+1)τ2_0 h_{i,j+1} + [j = i−1] iτ2_0(h_ii − conj h_ii)]/(2d−1).
> Exact counterexample to the printed form: k = 4, key (2,1) (actual h21 = −3√2/4; printed form is off by 3√2(1+i)/16).
> The code (step_39.py) uses the correct form, and with it Lemma D and the corollary hold for all k (CX4 review.md, WT2r result.md).

**PROVED — the Stabilization Lemma holds for all k, hence Conjecture 3 (weight-two part of M_k = (π²+6λ²)/24) holds for all k ≥ 2.**

Route used: **(c) induction on k with an explicit invariant of the tail tower** (top-order stabilization), powered by **(a) closed-form generating functions** for the initial data. Full write-up in `result.md`; exact data in `record_W2b.pkl` (all W, HigherPart, Gm for k ≤ 12).

## The proof

**Setup (WT2, re-verified).** The reduction is the unique ε-series solution of the triangular system on canonical keys (i,j), d=i+j:
W_ij = ε^d R_ij − (i+1)Gm_{i+1,j}DL1 − [i>j](j+1)Gm_{i,j+1}DL2 − [j=i−1]iF(Gm_iiDL1), Gm_ij = (1−eval_b)∘Int∘HigherPart(W_ij). DL1, DL2 have simple zeros at v=±i with τ1_0=(−1−i)/2, τ2_0=(−1+i)/2.

**Lemma A (pole-order bound).** HigherPart(W_ij) has ±i order ≤ 2d−2; Gm_ij ≤ 2d−3 (induction: Int lowers by 1, DL's zero lowers by 1, F preserves).

**Lemma B (top-order recursion).** h_ij^(k) := top (+i) coefficient of HigherPart(W_ij^(k)) obeys the *scalar* recursion h_ij = δ_{k,d}ρ_ij + [(i+1)τ1_0h_{i+1,j} + [i>j](j+1)τ2_0h_{i,j+1} + [adjacent: iτ2_0(h_ii − conj h_ii)]]/(2d−1) — the bound in Lemma A means only top coefficients of the feeds contribute.

**Lemma C (closed-form initial data).** ρ_pq = K_d(−1)^p(1+[p>q])binom(d,p), K_d = (−1)^d(2d−3)!!i^{d+1}/(2√2 d!) — derived from the master integrand (only the j=d term of the sqrt block contributes; uses (z′(i))^{2−2d}=2^{d−1}, U_d^{(m)}(τ)=i^d(−1)^m binom(d,m) via the finite-difference identity Σ_s(−1)^s binom(d,s)s^m(d−s)^{d−m}=d!(−1)^m, and the t↔1−t fold doubling off-diagonal keys).

**Lemma D (top-order stabilization).** (i) Off-diagonal: h_pq^(k)=ρ_pq for ALL k≥d. (ii) Diagonal: h_pp^(2p)=ρ_pp, h_pp^(k)=ĥ_p=(p+1)τ1_0ρ_{p+1,p}/(4p−1) for k≥2p+1. Proof by induction on level k (downward on d); constancy ⟺ four identities which follow exactly from Lemma C: (I1) K_{d+1}(d+1)=−i(2d−1)K_d; (I2) non-adjacent recursion (via binomial identities + τ2_0−τ1_0=i ⟹ (I1)); (I3) adjacent recursion (direct closed-form evaluation); (I4) Im(ĥ_p)=Im(ρ_pp) (via (I1)) so (I3) propagates past level 2p.

**Theorem (Stabilization Lemma, all k).** Degree-2 keys have order ≤ 2 (Lemma A), so their HigherPart IS the top coefficient: h_20^(k)=ρ_20=−√2i/4 ∀k≥2; h_11^(2)=ρ_11=√2i/4, h_11^(k)=ĥ_1=√2(1+i)/4 ∀k≥3. Applying (1−eval_b)∘Int gives exactly Gm_20(ε)=ε²/(1−ε)g20, Gm_11(ε)=−ε²g20+ε³/(1−ε)g11.

**Corollary (Conjecture 3, all k≥2).** Substituting into W_10 gives the exact ε-series identity (1−ε)(W_10−εR_10)=ε²Rw2 (c2=Rw2, c3=0, no other terms — verified symbolically), so [ε^k]SimplePart(W_10)=Rw2 ∀k≥2, and by the paper's reduction + Appendix dilog identity the weight-two part of M_k is ∫_b^a Rw2·log t dv = (π²+6log²2)/24. ∎

## Evidence (all exact, residual 0; "digits" = exact Q(√2,i) equality)

| check | range | result | step |
|---|---|---|---|
| (∗) == reducer W_ij | 238 items, k≤12 | MATCH, residual 0 | step_38 |
| Lemma A bounds | all data k≤12 | HOLD | step_38 |
| Lemma B recursion | all keys k=3..12 | MATCH | step_39 |
| ρ closed form == initial data | 49 keys, d=2..12 | MATCH | step_37/42 |
| R~ formula == buildintegrand | d≤6 all m | MATCH | step_41 |
| fold identity, U_d(τ), z′(i) | d≤6..8 | MATCH | step_41 |
| (I1)–(I4) | d≤16, p≤10, exact | HOLD | step_37/42 |
| Lemma D stabilization | all keys k≤12 | AS PREDICTED | step_42 |
| Gm_20=g20 k≥2, Gm_11 pattern | k=2..12 + k=13,14 | MATCH | step_40/43 |
| SimplePart(W_10^(k))=Rw2 | k=2..12 | MATCH | step_42 |
| (1−ε)Wcas=ε²Rw2 (symbolic ε) | all ε | TRUE | step_40/42 |
| kernel identity K_k=Rw2 | k≤15 (WT2) | MATCH | WT2 |

Key step files: step_31–33 (tail-tower data, polynomiality pattern), step_34 (ρ closed-form discovery), step_35–36 (F-branch/top recursion), step_37 (closed form + identities), step_38–39 (system + recursion + stabilization vs data), step_40 (Gm's + final identity), step_41 (initial-data formula lemmas), step_42 (full proof verification suite), step_43 (k=13,14 stabilization). Machinery: fastred.py/genintegrand.py/buildintegrand.py (from WT2).

**Paragraph for the paper:** We prove Conjecture 3 for all k. In the ε-generating form of the reduction, the IBP tower on the letters ±i closes on principal parts, and its leading coefficients at v=+i satisfy a scalar recursion whose initial data has the closed form ρ_pq=K_d(−1)^p(1+[p>q])binom(d,p), K_d=(−1)^d(2d−3)!!i^{d+1}/(2√2d!). Three algebraic identities — a recurrence for K_d and two binomial evaluations, one reflecting the v→1/v fold of diagonal keys — imply the top coefficients are level-independent: off-diagonal keys keep their initial value, diagonal keys stabilize one level later. Since the keys (2,0),(1,1) feeding the weight-two kernel have pole order ≤ 2, their higher parts are exactly these stable coefficients, yielding Gm_20=ε²/(1−ε)g20 and Gm_11=−ε²g20+ε³/(1−ε)g11. Substitution gives (1−ε)(W_10−εR_10)=ε²Rw2 identically in ε, hence the weight-two kernel of M_k equals Rw2 for every k≥2 and the weight-two part of M_k is (π²+6log²2)/24.