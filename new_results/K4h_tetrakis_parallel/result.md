# K4h complete — parallel h=1/768 / dps-270 production run for [d³]S of the tetrakis hexahedron is LAUNCHED and healthy

## 1. Parallelization (Task 1) — `k4h_driver.py`
`chamber_S`'s outer φ-grid (for model A a single interval [0, π/4] — `outer_splits(NA)` returns no internal splits) is built by the **same code path** as the original, cut into contiguous chunks of outer nodes; each worker runs a **verbatim copy** of `chamber_S`'s inner θ-quadrature block (`inner_sum`); chunk partials are summed in chunk order and ×48 at the end — the only difference vs single-process is cross-chunk summation order (roundoff level). `spawn` Pool, per-worker init builds model+node list once, atomic per-chunk JSON checkpoints (`mp.nstr(v, dps+guard+10)` — `str(mpf)` was tested NOT to round-trip), optional `guard` digits (production accumulates at 290 digits). Pitfall found & fixed: `__main__` guard is mandatory with spawn.

## 2. Validation (Task 2) — h=1/96, dps 80, WMIN 1e-45, LMIN 1e-35 (jobA_h384 convention), T=4.4
- Parallel 529 s vs sequential 2,819 s (5.3× on 6 workers); stats identical (zero_len 25, nbp 0, nsing 0; 608,006 inner evals).
- **Resume test: bit-identical** (deleted a checkpoint, reran, max diff 0.0).
- Agreement with single-process `chamber_S`, with a guarded 130-digit decomposition (step_04b) showing who owns the error:

| k | \|par−ref\| | E_par (parallel err) | E_ref (sequential err) | sig. digits agreed | ≤1e-75? |
|---|---|---|---|---|---|
| 0 | 2.0e-79 | 4.59e-80 | 1.54e-79 | 78.87 | ✓ |
| 1 | 1.0e-79 | 1.71e-80 | 8.29e-80 | 79.63 | ✓ |
| 2 | 7.0e-79 | 1.24e-78 | 1.94e-78 | 78.49 | ✓ |
| 3 | 5.82e-77 | 7.74e-77 | 1.36e-76 | 76.60 | ✓ |
| 4 | 2.47e-75 | 3.14e-75 | 5.62e-75 | 72.26 | ✗ (2.5×) |

k≤3 pass with 2–4 orders of margin. k=4 lands at 2.47e-75 vs the 1e-75 bound; the measured floor is the **sequential run's own error E_ref = 5.62e-75** (inner-sum cancellation on a 4.48e-3 value), so no parallel reordering at dps 80 can do better — the parallel error is smaller than the sequential one for every k, and the diff is fully accounted for. At production dps 270+20 guard this reordering class scales to ~1e-265, ~50 digits below the 1e-215 design floor.

## 3. Cost estimate (Task 3) — pilot at full production params, 6 workers
Range checked: weights reach WMIN=1e-215 at t*=5.752; **T=6.0** (w(6.0)=4e-276 ≪ WMIN) — the suggested 5.8–6.0 window confirmed. N=8,835 outer nodes (built at working dps; a dps-60 build undercounts 28% because tanh rounds to 1), 8,833 inner evals/node, per-node 75.0 s mean (loaded machine). **Single-process ≈ 184 h; 6 workers ≈ 30.7–31.6 h**; 369 chunks × 24 nodes ≈ 30 min each → ~5 min checkpoint cadence. Live: first 6 chunks took 1,804–1,811 s, log ETA 109,875 s ≈ 30.5 h — matches pilot.

## 4. Launch (Task 4) — DETACHED, running now
```
nohup nice -n 10 <project>/.venv/bin/python3 \
    <project>/tools/detach.py \
    <project>/.venv/bin/python3 \
    <project>/new_results/K4h_tetrakis_parallel/jobA_h768.py \
    > .../K4h_tetrakis_parallel/jobA_h768_out.txt 2>&1 &
```
(venv python used because PATH `python3` lacks mpmath.) Own session (probe-tested), effective nice **15** (`nice -n 10` on the session's 5, measured by the job via `os.nice(0)`). **Job PID 66142** (`jobA_h768.pid`); workers **66334, 66335, 66336, 66339, 66340, 66341, 66342**. Log: `bgA768.log`; checkpoints `ckpt_h768/chunk_XXXX.json` (6/369 done at write-up). On completion writes `modelA_h768.json` (same format as `modelA_h384.json`: Sdk ×5 at 270 digits, stats, HO, HI) + `modelA_h768_meta.json`, logging `[d^k]S` k≤4 and the d⁰ error vs 2√5/3. **Resume:** rerun the same command (or `python3 step_06.py`) — `run_job` skips existing checkpoints (validated bit-identical); `rm -rf ckpt_h768` for a fresh start.

## 5. Files
`result.md` (full report), `k4h_driver.py`, `jobA_h768.py`, `step_01.py`–`step_07.py` + `step_04b_guard.py` (step_02: splits/cutoffs/counts/benchmarks; step_03: three-way equivalence at h=1/24, max diff 6.5e-46 at dps 50; step_04: validation; step_04b: guarded decomposition; step_05: pilot/`pilot_h768.json`; step_06: launch/`launch_info.json`; step_07: monitor + result.md), logs `val96*.log`, `val96g*.log`, `bgA768.log`, results `val96_result.json`, checkpoints `ckpt_val96*/`, `ckpt_h768/`. Per instructions, the run was left running (not waited on); ETA ≈ 31 h wall under current machine load.