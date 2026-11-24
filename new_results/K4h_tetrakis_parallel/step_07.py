"""step_07.py -- K4h: collect monitor data and write result.md."""
import os, json, glob, subprocess, time
import mpmath as mp
mp.mp.dps = 100
CWD = "<project>/research/08_theory/K4h_tetrakis_parallel"
os.chdir(CWD)

log = open("bgA768.log").read().splitlines()
job_pid = open("jobA_h768.pid").read().strip()
workers = subprocess.run(["pgrep","-P",job_pid],capture_output=True,text=True).stdout.split()
nck = len(glob.glob("ckpt_h768/chunk_*.json"))
pilot = json.load(open("pilot_h768.json"))
val = json.load(open("val96_result.json"))
launch = json.load(open("launch_info.json"))
chunk_lines = [l for l in log if "chunk" in l and "done" in l]

# per-k digit agreement
Spar = [mp.mpf(s) for s in val["Spar"]]
Sref = [mp.mpf(s) for s in val["Sref"]]
Epar = [mp.mpf(x) for x in ["4.59e-80","1.71e-80","1.24e-78","7.74e-77","3.14e-75"]]
Eref = [mp.mpf(x) for x in ["1.54e-79","8.29e-80","1.94e-78","1.36e-76","5.62e-75"]]
rows = []
for k in range(5):
    d = abs(Spar[k]-Sref[k])
    sig = -mp.log10(d/abs(Spar[k]))
    rows.append((k, mp.nstr(d,3), mp.nstr(Epar[k],3), mp.nstr(Eref[k],3), mp.nstr(sig,4)))

eta_line = chunk_lines[-1] if chunk_lines else "(no chunks yet)"
md = f"""# K4h: parallel high-precision run for the tetrakis hexahedron's d^3 coefficient (production h=1/768, dps 270)

**Work dir:** `research/08_theory/K4h_tetrakis_parallel`. Inputs (read-only): `../K4e_tetrakis_d3` (k4c_driver.chamber_S, k4_engine, k4_quad, jobA_h384.py, modelA_h384.json) and `../K4g_tetrakis_z3/result.md` (why ~103 digits is not enough: PSLQ with 19 constants at 103 digits sees coefficients only to ~1e5.4, while the d^2 relation already has coefficients ~1.3e4).
Machine: Apple M4 Pro, 14 cores, 48 GB; python = the project venv (`<project>/.venv/bin/python3`, mpmath 1.3.0); the PATH `python3` (homebrew) has **no mpmath**, so the venv interpreter is used explicitly everywhere.
Labels: PROVED / NUMERICAL / CONJECTURE / PENDING. Every number below was printed by a step/job file here.

## Task 1 — parallel driver (`k4h_driver.py`)
`chamber_S` integrates over an outer phi-grid on [0, pi/4] (for model A `outer_splits(NA)` returns NO internal split points: a single interval, 6,907 nodes at dps 60 / **8,835 nodes at dps 270** — the count is dps-dependent because tanh-sinh u=tanh(...) rounds to exactly 1 at low dps and those tail nodes are then rejected by the `x<b` filter). The driver:
- builds the outer node list with the SAME code path as `chamber_S` (same splits, same `tanh_sinh_nodes`, same WMIN filter) — `build_outer_nodes`;
- cuts it into contiguous chunks of `chunk` outer nodes (production: 24 -> 369 chunks);
- each worker computes each node's inner theta-integral with `inner_sum`, a VERBATIM copy of the inner block of `k4c_driver.chamber_S` (breakpoint logic, LMIN skips, HI2 switch, g(th), weights w*w2 all identical);
- chunk partials are summed in chunk index order and multiplied by 48, exactly as `chamber_S` does at the end.
So the ONLY difference vs single-process is floating-point summation order across chunk boundaries (roundoff level). Multiprocessing: `spawn` context, `init_worker` builds the model + node list once per worker; checkpoint files `ckpt_*/chunk_XXXX.json` hold per-chunk partial sums as `mp.nstr(v, dps+guard+10)` strings (str(mpf) does NOT round-trip — tested; nstr(dps+10) does); writes are atomic (tmp + os.replace). Optional `guard` digits: workers accumulate at dps+guard (production uses guard=20 -> 290-digit accumulation; cost a few %, removes any reorder-noise concern: at dps 80 the measured reorder noise is ~3e-75 ~ 2.3e4*eps, scaling to ~1e-265 at dps 290, 50 digits below the 1e-215 design floor).
NOTE: all step scripts use the `if __name__ == "__main__":` guard — mandatory with spawn (an early unguarded version re-ran the whole script in every child and hung).

## Task 2 — validation at h=1/96, dps 80 (NUMERICAL)
Params (matching WMIN/LMIN per the jobA_h384 convention WMIN=1e-(dps-35), LMIN=1e-(dps-45)): **WMIN=1e-45, LMIN=1e-35**; TO=TI=**4.4** (step_02: for h=1/96 and WMIN=1e-45 the weights reach WMIN at t*=4.188, and w(4.4)=2.3e-56 << WMIN). chunk=32 -> 26 chunks, 805 outer nodes, 608,006 inner evals; 6 workers.
- Timing: parallel **529 s** vs single-process chamber_S **2,819 s** (5.3x on 6 workers, shared machine). Stats identical: zero_len=25, nbp=0, nsing=0.
- Resume test: deleted `ckpt_val96/chunk_0025.json`, reran `run_job` -> recomputed only that chunk (7 s); max |S_resume - S_par| = **0.0** (bit-identical).
- |S_par - S_ref| per k (bound 1e-75):

| k | par-ref | E_par = |guarded130 - par| | E_ref = |guarded130 - ref| | sig. digits agreed |
|---|---------|---------|---------|------------------|
"""
for (k, d, ep, er, sg) in rows:
    md += f"| {k} | {d} | {ep} | {er} | {sg} |\n"
md += f"""
- k <= 3 PASS with 2-4 orders of margin. k=4: |diff| = 2.47e-75, i.e. 2.5x over the stated 1e-75 bound. The guarded 130-digit decomposition (step_04b) shows the single-process chamber_S itself deviates from the exact sum of the same quadrature rule by E_ref = **5.62e-75** on k=4 (its value is only 4.48e-3; inner-sum cancellation, chunk-partial ratio ~23x, not tail cutoff: both schemes drop the identical ~1e-45 tail, and val96 vs K4e's modelA_h96/h384 differs exactly at that tail scale, 1.1e-44 for k=0). Since E_par (3.1e-75) < E_ref for every k, the parallel scheme is at least as accurate as the sequential one; |par-ref| is bounded below by ~E_ref minus correlated-error cancellation, so **no parallel reordering at dps 80 can reach 1e-75 on k=4** — the bound is unattainable there in principle; agreement is 72.3 significant digits on k=4 and >= 76.6 on the others (verdict recorded as measured in val96_result.json).
- Production impact: at dps 270 with guard=20 this reordering class scales with eps to ~1e-265, ~50 digits below the 1e-215 WMIN floor — irrelevant to the >= 200-digit target.

## Task 3 — production cost (pilot, step_05.py, NUMERICAL)
Production params: h=1/768, dps 270 (+20 guard), WMIN 1e-215, LMIN 1e-205, TO=TI=**6.0**. Range check (step_02): for h=1/768, WMIN=1e-215 the tanh-sinh weights reach WMIN at t*=5.752 (max over outer and inner interval lengths), and w(6.0)=4e-276 << WMIN — the suggested 5.8-6.0 window confirmed; T=6.0 loses no node with w > WMIN.
- Outer nodes N = **8,835** (built at working dps 270); inner evals per node: 8,833 typical (max 8,833; mean 8,508; a few near phi=0 are zero-length).
- Pilot: 60 strided nodes through the real 6-worker Pool: wall 772 s; per-node time min 0.0 / median 77.0 / mean 75.0 / max 94.8 s (machine shared with other jobs; idle benchmark gave 399 F_series-evals/s at dps 270).
- Estimate: single-process 662,765 s = **184 h**; 6 workers **30.7-31.6 h** (mean-based and throughput-based agree). Chunks: 369 x 24 nodes, ~30 min each per worker -> checkpoint cadence ~5 min wall.
- Live confirmation: first 6 production chunks took 1,803.9-1,811.0 s each (~30.2 min), log ETA after 6 chunks = 109,875 s = 30.5 h, matching the pilot.

## Task 4 — production launch (DETACHED)
Exact command (run from the work dir, via step_06.py):
```
nohup nice -n 10 <project>/.venv/bin/python3 \\
    <project>/tools/detach.py \\
    <project>/.venv/bin/python3 \\
    <project>/research/08_theory/K4h_tetrakis_parallel/jobA_h768.py \\
    > <project>/research/08_theory/K4h_tetrakis_parallel/jobA_h768_out.txt 2>&1 &
```
`detach.py` does os.setsid()+execvp: the job runs in its own session and survives this CLI session (probe-tested: child pid==launcher pid, own sid). Niceness: launched `nice -n 10` on top of this session's 5 -> effective nice **15** (measured by the job itself via os.nice(0), first log line); the 6 spawn workers inherit it.
- **PIDs**: job (session leader) **{job_pid}** (also in `jobA_h768.pid`); worker children: **{', '.join(workers)}** (6 Pool workers + 1 spawn resource tracker).
- **Logs**: `bgA768.log` (progress, one line per chunk with running ETA; stdout mirror `jobA_h768_out.txt`).
- **Checkpoints**: `ckpt_h768/chunk_XXXX.json` (idx, n, acc[5] at 300 digits, stats, elapsed), atomic.
- **Outputs** (on completion): `modelA_h768.json`, SAME format as `modelA_h384.json` ({{"Sdk": [5 strings at 270 digits], "stats": ..., "HO": ..., "HI": ...}}), plus `modelA_h768_meta.json` (dps, T, WMIN, LMIN, guard, workers, chunk, N, elapsed, d0 error vs 2sqrt5/3).
- **Time estimate**: ~31 h wall under current machine load (pilot 30.7-31.6 h; first-round live ETA 30.5 h), i.e. completion ~2026-10-28 23:30 local if the load persists; faster if the machine frees up.
- **How to resume**: just rerun the exact command above (e.g. `python3 step_06.py`, or the nohup line). `run_job` scans `ckpt_h768/`, skips every chunk with a valid JSON, recomputes the rest (corrupt files are ignored and recomputed). Validated: resume after deleting a chunk reproduced the final sums bit-identically (max diff 0.0, step_04). To restart from scratch instead: `rm -rf ckpt_h768` first.
- **Monitor**: `tail -f bgA768.log`; completion marker `jobA768 DONE`; the log also prints [d^k]S for k<=4 at 250 digits and the d^0 error vs 2sqrt5/3 (expected ~1e-215 or better).

## Status at write-up (step_07, {time.strftime("%Y-%m-%d %H:%M:%S")})
- {nck}/369 chunks checkpointed. Last log line: `{eta_line}`
- Job alive: {'yes' if workers else 'CHECK'} (pid {job_pid}, workers {len(workers)} child processes).

## Step files
- step_01.py recon of inputs; step_02.py outer splits, tanh-sinh weight cutoffs (t* bisection), node counts, h384 calibration, F_series rate benchmarks;
- step_03.py cheap three-way equivalence h=1/24 dps 50 (orig chamber_S vs sequential node loop vs parallel; max diff 6.5e-46, roundoff level) + str/nstr round-trip and node-list determinism tests;
- step_04.py validation h=1/96 dps 80 (parallel + resume test + sequential reference; logs val96.log/val96_out.txt, val96_result.json);
- step_04b_guard.py guarded 130-digit error decomposition (E_par vs E_ref; ckpt_val96_g130/);
- step_05.py production pilot (pilot_h768.json); step_06.py production launch (launch_info.json); step_07.py this summary;
- k4h_driver.py (parallel driver), jobA_h768.py (production script), detach probe (removed), worker_pids.txt.
"""
open("result.md","w").write(md)
print(md[:3000])
print("...")
print("result.md written,", len(md), "chars")
print("chunks done:", nck, "| job pid:", job_pid, "| workers:", workers)
