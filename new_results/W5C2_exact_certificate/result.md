W5C2 complete — EXACT certificate for W5 = N5 + R5_min3, all 8 systems verified exactly (second fresh verify run also passed).

# W5C2 — exact certificate for M4's weight-5 part

## Verdict: **EXACT / COMPUTER-ASSISTED** — all 8 systems (4 Galois sectors × {D_rat, D_sq}) verified exactly; nothing partial.

For D = W5 − (N5 + R5_min3) as exact level-8 H-word vectors (46656 admissible weight-5 words), explicit rational vectors y were produced with **Σᵢ yᵢ Rᵢ = Φ_s(D) EXACTLY (Fractions, residual 0)** in each of the four Galois sectors for both components — the exact analog of A3h5's W4 certificate at weight 5.

## Per-system evidence (all printed by my code; logs in progress.log)

| sector | rows×cols | rank = \|S\| = \|C\| | S families | Dixon steps | y height num/den (digits) rat | sq | p2 checks | EXACT residual (rat, sq) |
|---|---|---|---|---|---|---|---|---|
| (+,+) | 83147×12218 | 12156 | engine 12120 + framew5 36 | 577 | 2078 / 2065 | 2081 / 2069 | ✓ | 0, 0 |
| (+,−) | 76044×11164 | 11101 | engine 11065 + framew5 36 | 563 | 2023 / 2012 | 2027 / 2018 | ✓ | 0, 0 |
| (−,+) | 76181×11164 | 11041 | engine 11041 (all) | 553 | 1991 / 1981 | 1991 / 1982 | ✓ | 0, 0 |
| (−,−) | 82557×12110 | 11999 | engine 11999 (all) | 552 | 1990 / 1979 | 1989 / 1980 | ✓ | 0, 0 |

- Second-prime (p2=16777199) confirmations: rank(S)=|S|, det(R_{S,C})≠0, D∈span(S) — PASSED all 4 sectors.
- Rank matches W5C1's Wiedemann lower bounds (pp 12156 = minpoly degree).
- Cross-checks A·Y = B mod p2 passed for all 4 solves.

## Proof chain (each step verified by code)
1. **Trimming** (LeanEchelon sparse GE mod p1): S independent rows, C pivot columns, D reduces to 0 (membership).
2. **Exact solve**: custom Dixon p-adic lifting (flint nmod LU = one n³ inverse + lifting, early-terminating rational reconstruction with stability repeat), y·R_{S,C} = D_C over Q, both RHS at once.
3. **Exact verification** (Fractions): Σ yᵢRᵢ = Φ_s(D) on ALL sector columns — residual 0, all 8 systems. Every used sector row = (1/λ)·fold of an unfolded source row (recomputed from rels_*.pkl / validated_w5.pkl), λ matches.
4. **Fold/unfold lemma**: u = Σ_s (1/4)Σ_g χ_s(g) g(Φ̃_s(u)) verified on all 46656 unit vectors; Σ_s n_s = 46656 ✓.
5. **Exact unfolded identity** (step_05, Fractions): D = Σ_s Σᵢ (yᵢλᵢ/4)Σ_g χ_s(g)(g·rᵢ) on ALL 46656 words — **residual 0**, both components.
6. **Value argument**: (g·r)·V = r·(gV) = g(r·V) = 0 (r proved relation with rational coeffs; g ∈ Gal(Q(ζ8)/Q) acts on cyclotomic periods) ⟹ D·V = 0 ⟹ **W5 = N5 + R5_min3**. No Galois set-closure needed (step_06: 36016/40392 transforms are literally in the files anyway).
7. **Standalone verify_W5.py** (fresh process, reads only certificates + W5C1 sources): orbit structure OK; lemma OK; all 8 sector equations residual 0; unfolded assembly over all 4 sectors mod 2^61−1 residual 0 (rat 25406, sq 28535 words). VERDICT: EXACT. Run twice, exit 0, ~60 s.

## Row families used (for the planner's provenance confirmation)
- **engine/standard (std)**: regularised double shuffle (Ihara–Kaneko–Zagier) + regularised distribution d∈{2,4,8} (Zhao) — pp 9835, pm 9019, mp 8980, mm 9681 rows.
- **engine lifts**: shuffle lifts of proved lower-weight relations — lift1 (2248/2007/2037/2275), lift2old (20/8/12/22), lift3old (4/6/6/8), lift4old (6/11/5/10).
- **frame-w5** (288 weight-5 frame double-shuffle rows, W5C1-validated to 2.8e-113): **36 used in (+,+), 36 in (+,−), 0 in (−,+)/(−,−)** — consistent with W5C1's certified necessity of frame-w5 in (+,+)/(+,−); engine alone spans D in the (−,·) sectors.
- **frame-234 lifts**: NOT used in any sector.
Rows are identified in the certificates by (family, index-in-file, label) for regeneration from rels_std/lift1/lift2old/lift3old/lift4old.pkl and validated_w5.pkl.

## Timings (code-printed)
step_01b trim: 679/992/385/371 s; step_02: 249/251/232/229 s; step_03 Dixon: 963/805/776/905 s (nmod inv(12156)=406 s); step_04 verify+cert: ~20 s each; step_04a 34 s; step_05 192 s; verify 60 s. Pipeline compute ≈ 2.3 h (+ development/benchmarking; LeanEchelon written after sparse-heap/dense/RREF prototypes proved too slow).

## Files (work dir `new_results/W5C2_exact_certificate`)
- `certificate_{pp,pm,mp,mm}.json` (~90 MB each): rank, |S|, primes, Dixon stats, heights, family counts, per used row {S_index, fam, famidx, tuple, y_rat, y_sq, src{ufam, uidx, label, λ}}.
- `verify_W5.py` (standalone), `result.md`, `summary.md`, `progress.log`.
- Steps: step_01b.py (trim), step_02.py (p2 checks), step_03.py (Dixon), step_04.py (exact verify + certificates), step_04a.py (fold replay/preimages), step_05.py (exact unfolded identity), step_06.py (Galois set-closure info); modules w5c2_lib.py, w5c2_fold.py, w5c2_dixon.py; unit tests step_0N_unittest.py (all echelons validated vs flint).
- `work/`: step01/02_*.pkl (S, C, checks), yexact_*.pkl, feed_*.pkl, foldmap_*.pkl.

## Caveats (honest)
- Depends on W5C1's D being correct (validated numerically to 1e-110 there) and on the row families being proved relations (planner's premise: engine families are the cited theorems; frame-w5 are W5C1's validated rows).
- The Galois action on values is the standard conjugation of (regularised) cyclotomic iterated integrals — the same action W5C1's sector construction is built on.
- One bug was caught and fixed during the run (stale rows_*.pkl vs feed_*.pkl ordering → singular A); the corrected pipeline then verified exactly end-to-end.
## ERRATUM (planner, 2026-11-22)

Step 6 above (the "value argument" for u = 3, 5) is false as stated: the image of a relation under √2 ↦ −√2 need not be a relation
(weight-1 counterexample and the affected rows are in `../W5G1_frame_galois/result.md`). The exact vector identity of this certificate
is correct, but as written it does not prove W₅ = N₅ + R₅. See W5G1 for the repair and for what remains open (the images of 37 frame
rows). This erratum is from an audit that is itself being refereed (`../W5G1R_referee_frame_galois/`).
