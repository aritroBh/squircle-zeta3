EXACT CERTIFICATE PASSED — Conjecture 4 PROVED (c0 = c1 = 0, hence W4 = R).

## What was accomplished

The resumed run executed the planner's fast Stage-2 plan and completed the full certificate:

**Certificate.** `D_X = Σ_i y_i^(X)·T0_{S_i} + Σ_j z_j^(X)·F_j` for X ∈ {rat, sq}, exact rational identities on the 5184 admissible weight-4 level-8 words, with 5040 T0 rows + 28 frame lifts = 5068 rows total (4810 std, 223 lift1, 7 lift2old, 28 frame). No Kummer rows, no R1–R5, no lift3old, no P4/π⁴ term.

**Pipeline (all timings printed by my code):**
- Stage 1 (planner, prior run): z for the 28 frame lifts mod 8 primes, CRT×7 + hold-out agreement; heights tiny (denominators divide 684 / 1824, max |num| 1729451).
- `step_29`: 5042 independent T0 rows S via transpose-RREF (5184×17625) mod p1, 129 s.
- `step_30`: G_X = D_X − Σ z_j F_j exactly (supports 3922/4016; den lcms 1368, 36480).
- `step_31`: B = S[:,pp]^T nonsingular mod p1 (13.6 s).
- `step_32`: exact `fmpz_mat.solve` of the 5042×5042 system, both RHS at once, 74.8 s; y heights max 381 digits.
- `step_33`: **exact Fraction verification y^T·S = G_X on ALL 5184 columns — residual 0, both components.**
- `step_34`: `D_certificate.json` (7.2 MB).

**Standalone verifier** (`verify_D.py`, copied from A3h4 with only three path adaptations: two `sys.path.insert`, the certificate path, and the log destination — `A3h4_frame_ds/progress.log` was not writable): it rebuilds *every* row from first principles (`make_relations`, lower bases asserted ⊆ generated relations, frame rows re-derived from provenance via candidates/convert/identity with domain assertions, target rebuilt from cmzv source) — reads none of my pickles. Clean fresh run, exit 0, 6.9 s:
```
EXACT rebuilt D component=0 support=3117 residual_nonzero=0 c0=0
EXACT rebuilt D component=1 support=3480 residual_nonzero=0 c1=0
VERDICT: D=(0+sqrt(2)*0)*pi^4; Conjecture 4 PROVED
```

**Sanity numerics (step 4, not part of the proof):** A3f4's validated evaluator (roots recomputed after raising prec to 900 bits/cap 2600, N=2400, ~270-digit accuracy; two-N check 3e-232): 10 random raw weight-4 frame rows and 10 random certificate lifts vanish to ~1e-269..1e-271; both full stage-1 combinations (supports 3275/3420 words) vanish to 2e-267 / 7e-266.

**Provenance (step 5):** regularised double shuffle (Ihara–Kaneko–Zagier), Zhao's regularised distribution d∈{2,4,8}, weight-1 branch identities, shuffle lifts of proved relations, and function-level double shuffle in real Möbius frames with enforced absolute-convergence domains (no branches, no continuation) per PLANNER_AUDIT.md. Caveats flagged in result.md (frame code path audited not re-implemented; std label indexing immaterial to soundness).

All artifacts in `.../08_theory/A3h5_exact_certificate/`: `result.md` (verdict first), `progress.log`, `D_certificate.json`, `verify_D.py`, `step_29.py`–`step_36.py`, stage pickles, `verifier.log`, `step_36.out`, with rerun commands listed in result.md.
## Planner verification (W4 side session, 2026-10-22 15:3x)
1. The verifier copy differs from A3h4_frame_ds/verify_D.py only in sys.path lines, the certificate path and the log
   destination. Checked with `diff`.
2. Fresh planner re-run of `verify_D.py` (`planner_rerun_verify.out`): residual 0 in both components, c0 = 0, c1 = 0,
   "Conjecture 4 PROVED", exit 0. Rows used: std 4810, lift1 223, lift2old 7, frame 28 (10 of weight 2, 9 of
   weight 3, 9 of weight 4).
3. The target is the paper's statement:
   - `w4_symbol.conjecture_terms()` (used by `rebuild_target`) was compared term by term with
     \begin{conjecture}\label{conj:W4} in paper/squircle_note/main.tex (lines 121–136). All 33 coefficients match,
     with L = log(1+√2) and λ = log 2 (main.tex line 53).
   - The paper's formula, evaluated directly with mpmath, equals W4 = 0.79578505983224081132… to the 48 digits
     available (|R − W4| = 2e-50, the truncation of the W4 literal).
4. Stage 1 (z) was reconstructed independently by the planner (`planner_probe_z.py`): CRT over 7 primes, with an 8th
   as hold-out.
5. Frame-row provenance: `../A3h4_frame_ds/PLANNER_AUDIT.md` (absolute-convergence domains, no branch choices, path
   convention), plus independent numerics (weight 3 to 1e-139 against A3f6 values).

STATUS: Conjecture 4 is PROVED, c0 = c1 = 0, W4 = R exactly. The proof uses the following inputs:
- (i) the proved G-expression of W4 (A3_W4_proof);
- (ii) theorems: regularised double shuffle and distribution for level-8 multiple polylog values, shuffle products of
  convergent iterated integrals, and function-level double shuffle of absolutely convergent multiple polylog series;
- (iii) correctness of the exact conversion code (cmzv path composition and Cayley pullback; frame_ds). It is audited
  and validated numerically (conversions to 1e-50, every used row to ≤ 1e-266), but not formally verified.
