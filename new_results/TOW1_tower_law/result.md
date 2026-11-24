# TOW1 — the tower law at v = i: weight-3 theorem for all k, polynomial-in-k law at every weight

## VERDICT

**A: PROVED.** The tower law at depths δ ≤ 2 is proved for **all degrees d and all k** (k-ranges
explicit: off-diagonal keys from k = d, diagonal keys from k = 2p+1). Every step is either proved
in writing (proof_tower.md) or verified by exact code (verify_tower.py, 12/12 PASS in 93 s).
Together with WT3's PROVED reduction theorem this proves the **weight-3 theorem W3(k) for all
k ≥ 3** (kernels R_20(k;v), R_11(k;v) and W3(k) = W3' − (k−3)(ζ(3)/8 + π²log2/24 + log³2/12)).

**B (all depths): REDUCED, with the mechanism identified.** The Newton-coefficient fixed-point
formulation (proved) reduces the full tower law to an explicit infinite triangular family of
identities; the δ ≤ 2 instance is Task A. The pattern of the top Newton coefficients
(e(2,1) = −e(0,0)/(2d−3); depth-4 top slope under investigation) and the all-depth initial-data
problem (Taylor coefficients of one analytic factor A(v) at v = i) are identified; see §7 of
proof_tower.md. Exact evidence: the tower law holds in every one of the 64 (n,δ) cells tested
(WTS-theory) and all our checks.

## Task A: what is proved (proof_tower.md has everything)

| item | status |
|---|---|
| (T1) depth-form recursion, all keys/depths | PROVED in writing + exact check (2207 cells, k ≤ 13) |
| Initial data ρ(δ), δ = 0,1,2, all keys | PROVED (S(t*), S'(t*) single/double-term evaluations; analytic-factor Taylor; Corollary S3 rational (S'/S)(t*)) + check (98 cells, d ≤ 10) |
| Fold identity: child diagonal feed = −p·DL2·[Gm(v)−Gm(1/v)] | PROVED (DL1(1/v)/v² = −DL2(v)) + check |
| Tower law δ ≤ 2 ⟺ identities E0,E1,E2a,E2b,E3 | PROVED (Newton-coefficient induction) |
| E0 (δ=0 fixed point) | PROVED (WT2b/CX4/CX5) |
| E1 (δ=1 identity) | PROVED all d (non-adjacent: hand; adjacent: exact rational simplification to 0, prove_E1_adjacent.py) |
| E2a (δ=2 slope eq) | PROVED all d: with the closed-form slope e(2,1) = −e(0,0)/(2d−3) it IS E0 |
| E2b (δ=2 value identity) | PROVED all d (exact rational simplification to 0, prove_E2b_nonadjacent.py, prove_E2b_adjacent.py) |
| E3 (diagonal blindness, δ ≤ 2) | PROVED all p (prove_E3.py; three identities → 0) |
| The tower law δ ≤ 2 itself | the above + exact check: all keys, k ≤ 14 (0 bad cells) |
| W3(k) for all k ≥ 3 | PROVED: tower law (deg-3 keys) ⇒ affine Gm law (slopes verified exactly, k ≤ 14) + WT3's reduction theorem |

## Task B: status
REDUCED to the infinite Newton-identity family (proved formulation); δ ≤ 2 solved (Task A);
mechanism and next targets identified (proof_tower.md §7). Not claimed beyond δ = 2.

## Files
- proof_tower.md — full proof document
- verify_tower.py — 12/12 PASS in 93 s (reruns every identity: data checks V1–V8, all-d symbolic proofs V9)
- prove_E1_adjacent.py, prove_E2b_nonadjacent.py, prove_E2b_adjacent.py, prove_E3.py, proof_common.py — the all-d symbolic proofs
- towlib.py — library; step_*.py — harness log of all computations
