result = r'''# TOW1 — the tower law at v = i: weight-3 theorem for all k, polynomial-in-k law at every weight

## VERDICT (live, updated as the job progresses)

**A (tower law at depths δ ≤ 2, all degrees): PROVED modulo three explicit identities.**
The proof is COMPLETE in structure: the tower law at depths δ ≤ 2 for all degrees d and all k
is reduced, by a proved recursion and proved closed forms, to three families of exact identities
(E1, E2b, E3 below). Each identity is **VERIFIED BY EXACT CODE for d ≤ 11** (k ≤ 14 for the
tower law itself). Upgrading the three identities from VERIFIED to PROVED (they are explicit
rational identities in (d,p)) is the remaining work; the slope family is closed-form
(e(2,1) = −ρ(0)/(2d−3)) and the δ=2 *slope* equation is shown to be **equivalent to WT2b's
already-proved δ=0 fixed point**, so only E1, E2b, E3 remain.

**B (all depths): REDUCED** — the Newton-coefficient fixed-point formulation (this job, proved)
identifies exactly the infinite family of identities needed; pattern for the slope at depth 4
also found (see proof_tower.md §7). Not yet attempted beyond δ ≤ 2.

## The theorem and the corollaries

**Tower law (Conjecture T of WTS-theory), depths δ ≤ 2.** For every key (p,q) (p ≥ q),
d = p+q, let h_{pq}^{(k)}[r] be the coefficient of (v−i)^{−r} in HigherPart(W_{pq}^{(k)})
(r = 2d−2−δ). Then, with the k-ranges made explicit:
- off-diagonal keys (p > q): h_{pq}^{(k)}[2d−2−δ] = Σ_{s=0}^{⌊δ/2⌋} binom(k−d, s) e_{pq}(δ,s)
  for ALL k ≥ d (the initial level k = d included);
- diagonal keys (p,p): h_{pp}^{(k)}[2d−2−δ] = Σ_s binom(k−2p−1, s) ê_{pp}(δ,s) for k ≥ 2p+1
  (the initial value at k = 2p is exceptional, WT2b Lemma D(ii));
with e(δ,0) = ρ(δ) (the initial data, closed forms proved), e(2,1) = −e(0,0)/(2d−3),
ê(δ,s) = (p+1)Σ_{δ'} τ1_{δ−δ'} e_{p+1,p}(δ',s)/(4p−δ'−1).
Degrees: constant at δ = 0, 1; affine at δ = 2 — i.e. degree exactly floor(δ/2).

**Corollary (WT3, the weight-3 theorem).** The affine Gm law of WT3 follows: the degree-3 tower
(keys (3,0), (2,1)) at depths ≤ 2 gives Gm_30^(k), Gm_21^(k) affine in k for k ≥ 3 with slopes
s30, s21 as in WT3; by WT3's **reduction theorem (PROVED for all k)**, W3(k) = W3' − (k−3)
(ζ(3)/8 + π²log2/24 + log³2/12) for all k ≥ 3, with the explicit kernels R_20(k;v), R_11(k;v).
(Verification of the s30/s21 values: exact code, see verify_tower.py.)

**Corollary (WTS kernel law at weight 3).** By WTS-theory's Lemmas 2–3 (proved there), the
tower law at depth ≤ 2 for all degrees also gives the weight-3 case of the kernel law.

## What is PROVED vs VERIFIED (details in proof_tower.md)

| item | status |
|---|---|
| (T1) depth-form order-m recursion, all keys/depths | PROVED (projection of WT2/WT2b fixed-point system P1); exact-checked 2207 cells k ≤ 13 |
| Initial data ρ(δ), δ = 0,1,2, closed forms | PROVED (Lemma C + S(t*), S'(t*) single/double-term evaluations, analytic-factor Taylor); exact-checked 98 cells d ≤ 10 |
| e(2,1) = −ρ(0)/(2d−3) satisfies slope eq ⟺ δ=0 fixed point | PROVED (algebra, this job) |
| Tower law δ ≤ 2 ⟺ identities E1, E2b, E3 (+ WT2b's E0) | PROVED (Newton-coefficient induction, this job) |
| E1 (δ=1 identity), E2b (δ=2 value identity), E3 (diagonal blindness δ ≤ 2) | **VERIFIED d ≤ 11** (exact); proof for all d in progress |
| Tower law δ ≤ 2 itself, all keys | VERIFIED k ≤ 14, d ≤ 10 (0 bad cells) |
| W3(k) for all k | follows from A + WT3 reduction theorem (WT3: PROVED) |

## Files
- proof_tower.md — full proof document
- verify_tower.py — reruns every exact identity (PASS/FAIL per lemma)
- towlib.py — library (tower data, recursion, e-families, closed forms)
- step_*.py — the harness log of every computation
'''
open('result.md','w').write(result)
print('result.md written', len(result))
