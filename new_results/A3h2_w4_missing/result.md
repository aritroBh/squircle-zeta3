# A3h2 — the 29 missing weight‑4 level‑8 relations (final report, resumed run)

Workdir: `new_results/A3h2_w4_missing`. **Status: steps 1–3 complete and verified; step 4 skipped by planner addendum (it is A3h5's job; A3h4 already has the mod‑p version); step 5 only partially executed — the run was interrupted a second time before the membership check could run. The incomplete piece is clearly marked below; nothing unverified is claimed.**

## Verdict
- The 29 missing relations **exist numerically and are certified**: found by single‑embedding (k=1) LLL against 258‑digit certified values of the 110 free words of the T3 quotient; every candidate re‑verified at 2× precision to **1e‑679** (N=3800, prec 2500 bits); the 29 relation vectors are mod‑p independent (rank 29); gap to accidental dependencies ≈ **1e46**. Per the rules, these remain **numerical hints, not proofs**, until the exact ℚ(√2) certificate (A3h5) closes them.
- The brief's "Galois trick" is **only half true** (step 1): it is valid for form‑level engine rows at all k, but Rogers/R1–R5 hold only at k=1 and k=7; hence the full search had to be single‑embedding.

## Step 1 — Galois trick validation (printed; prec 700 bits, N=800)
σ_k: letter j↦k·j mod 8, −1 fixed; on coefficients √2↦(+,−,−,+)·√2 for k=(1,3,5,7).

| test | k=1 | k=3 | k=5 | k=7 |
|---|---|---|---|---|
| Rogers five‑term (36 wds) | 1.4e‑120 | **11.08** | **11.08** | 1.4e‑120 |
| R1 | 8.8e‑119 | **22.8** | **22.8** | 8.8e‑119 |
| R2 | 3.2e‑118 | **185** | **185** | 3.2e‑118 |
| R3 | 3.3e‑118 | **178** | **178** | 3.3e‑118 |
| R4 | 5.4e‑118 | **42.2** | **42.2** | 5.4e‑118 |
| R5 | 8.8e‑118 | **250** | **250** | 8.8e‑118 |
| std T3 rows (w4) | 5.0e‑120 / 1.1e‑120 | 5.0e‑120 / 2.0e‑121 | same as k=3 | = k=1 |
| shuffle identities | ≤1e‑188 | ≤1e‑188 | ≤1e‑188 | ≤1e‑188 |
| random ℤ‑combo (control) | 280.6 | 105.6 | 105.6 | 280.6 |

Identified residual: five‑term at k=3 fails by exactly **−4·π·i·log(1+√2)** (ratio to 40 digits). Conclusion printed in `progress.log`: trick holds for engine/form‑level rows, **fails for Rogers/R1–R5 at k=3,5** (branch shifts under √2↦−√2); k=7 (complex conjugation) is fine.

## Step 2 — values (printed)
- `vals110.pickle`: 110 free words × k∈{1,3,5,7}, acb prec 1200/cap 2400, N=2000. Convergence certificates: max|N2000−N1000| = 3.6e‑190 (k=1,7), 1.5e‑192 (k=3,5); max|N2000−N1400| = 2.2e‑264 ⇒ **certified ≥258 digits** (self‑consistency ~369 d). Worst words (7,7,7,0),(7,7,5,0).
- Independent nested `mp.quad` for G((7,4,7,7);1): matches engine to **1.62e‑38** (596 s).
- Evaluator cross‑validated vs A3f4's certified `vals_hi.pickle`: 5.6e‑151 on the 32 weight‑3 free words.

## Step 3 — the 29 relations (printed; DONE, not redone)
Single‑embedding LLL, C=1e300, on the k=1 values:
- exactly **29** kept relations (residual ≤ threshold at re‑check; gap 1e46 to accidents);
- each re‑verified at 2× precision: residual **1e‑679** (N=3800, prec 2500);
- mod‑p rank of the candidate set = **29**;
- heights (max |coeff|) from **1 to 264**; support sizes **6 to 109** words;
- only relation N0 is form‑level (invariant at all four k); the other 28 are "exotic" (fail k=3,5 like Rogers).
- N0 (printed from `w4_missing_relations.json`): the 6‑term, depth‑4, letters⊆{3,6,7} relation
  `G(7,6,3,3) − G(7,6,3,6) − G(7,6,6,3) + G(7,6,6,6) + G(7,7,6,3) − G(7,7,6,6) = 0` (all coeffs ±1).

Artifacts: `w4_missing_relations.json` (word→int coeff, all 29), `cands29.pickle`, `cands29_reduced.pickle`, `vals110_k1_700.pickle`.

## Step 4 — D target: SKIPPED per resume addendum
A3h4 (`A3h4_frame_ds/result.md`, Codex) showed mod both primes 1000000009/1000000033 that D lies in span(T3 + proved frame double‑shuffle lifts) + ℚ·P4 (T3 rank 5074→5082, quotient 110→**102**; also a T0‑based rational‑only route, quotient → 28+1). The exact ℚ(√2) certificate `verify_D_target.py` is **A3h5's job** — not claimed here.

## Step 5 — structure: PARTIAL (interrupted)
Done and printed: N0's depth/letter profile (above); support/height ranges for all 29 (6–109 words, heights 1–264); the full per‑relation depth/letter table is derivable offline from `w4_missing_relations.json` in one line (`depth = #{letters ≠ −1}`).
**Not executed** (second interruption hit while loading `validated_w*.pickle`, which needs `sys.path` to include `paper/squircle_note/scripts` for the `qsqrt2` module — that missing import was the last error): the membership test of each relation in span(T3 + A3h4 frame lifts) mod p1, and the dihedral/Galois orbit reduction.
Ready‑to‑run recipe (validated components all exist): load `rref_T3_1000000009.pkl` (A3h) for the 110‑dim nf map; nf each of the 29 relations; lift frame rows exactly as in `A3h4_frame_ds/planner_check_membership.py` (`cmzv.mul` + `relations.primitive` over all complementary words); test membership with `nmod_mat.rank`. **Expected outcome (A3h4's printed ranks, their numbers not mine): the frame span inside the 110‑dim quotient is 8‑dimensional, so an 8‑dim subspace of the 29‑dim missing space is absorbed by span(T3+frames) and ≈21 of the 29 basis relations should test "outside".** This must be confirmed by running the recipe; it is not a result of this run.

## Files
`evsetup.py` (trap‑safe evaluator setup), `rref_driver.py` + `rref_T3_{p}_b.pkl` (full T3 RREFs incl. nf(P4,Drat,Dsq)), `vals110.pickle`, `vals110_k1_700.pickle`, `cands29*.pickle`, `w4_missing_relations.json`, `step_01..57.py`, `progress.log`, `result.md` (in‑progress draft; replace with this report).

## Rerun commands
```
PY='<project>/.venv/bin/python'
cd new_results/A3h2_w4_missing
$PY rref_driver.py 1000000009 b     # rebuild T3 rref artifacts
# step-5 membership (one script): sys.path += [A3c_claude_elementary/Q, paper/squircle_note/scripts];
#   nf via rref_T3_1000000009.pkl; frames via A3h4_frame_ds/planner_check_membership.py pattern;
#   rank tests with flint nmod_mat.
```
Caveat repeated: all 29 relations are LLL‑derived numerical identities (certified to 680 digits, mod‑p rank 29); they are proof targets, and the exact D‑certificate is A3h5's task.
## Planner check of step 5 (W4 side session, 2026-10-22)
`planner_check_29_in_frames.py`, run with A3h's T3 RREF at p = 1000000009 and p = 1000000033, gives identical results
at both primes:
- In the T3 quotient (110 dimensions), the 29 relations have rank 29, and T3 + the 29 has quotient exactly 81 = 3⁴.
- A3h4's proved frame lifts span 8 dimensions there. Their span lies entirely inside the span of the 29
  (rank of frames + 29 = 29).
- No single basis relation of the 29 lies in the frame span on its own. This depends on the basis chosen.
