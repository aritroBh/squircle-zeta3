# W5G3: certificate store for the 74 conjugated frame relations -- DONE

Store: `W5G3_cert_store/store/` (148 files: `<name>.core.jsonl.xz`, `<name>.lit.jsonl.xz` for 74 targets), `MANIFEST.txt` (size + sha256 of every file).
Size: **14.097 GB total = core 3.103 GB (a_j, W, t_k) + lit 10.994 GB (literal-row coefficients)**; xz preset 3. Limit was 25 GB.
Format: `FORMAT.md`. Checkers: `check_store_example.py` (smoke test only), `compare_aj.py`.

## Results
- t03 (4 sectors, exact Dixon), rerun 09:33-10:02: same Dixon steps and digit sizes as first run.
- t04 (exact weight-5 identity on all 46656 words): **0 mismatches for 74/74** (45 "mismatches = 0" lines in `t04.run1.log`, of which 44 are
  kept; 30 in `t04.log`, which includes the redone u3.f68; 44 + 30 = 74). `t04.log` ends "ALL 74 targets ... verified".
- t05 (exact weight-4 identity W = sum t_k G_k on 5184 words): **0 mismatches for 74/74** (`t05.log`, "ALL W in the proved weight-4 span (exact): True", 6869 s).
  `w4cert_summary.pkl` equals the first run's (n_gens, frame lifts, frame rows, max digits, identical dicts).
- a_j: **identical to W5G2 `cert_summary.pkl` for 74/74 targets**, compared as exact Fractions keyed by frame-row identity (`compare_aj.out`).
  W support and frame-term counts also equal first run for all 74.
- Smoke test `check_store_example.py` (written from FORMAT.md only, reads the store files and the on-disk row files, no W5G3/W5G2 scripts):
  u3.f4, u3.f68, u5.f113: I5 and I4 mismatches = 0 (`smoke_*.out`). u3.f4 is a run-1 target, u3.f68 was redone in run 2.
- `xz -t` on all 148 files: no errors (`xz_test.out`).

## Interrupted run and resume (machine reboot ~12:13)
- Run 1 (t04, 10:15 -> killed ~12:13) verified and wrote targets 1-45 (u3.f4 ... u3.f68 in `t04.run1.log`). u3.f68's lit file was 0 bytes.
- Run 2 (`RESUME=1`, 12:18 -> 13:33, 4500 s): for each target both files must pass certio.verify_pair (complete xz decode, every JSON line parses,
  meta/frame/W/rows present, term counts equal meta counts; sha256 of each decoded stream logged). 44 targets passed and were skipped; u3.f68 and the 29 never-reached
  targets were recomputed and verified exactly.
- **Verified in run 1** (exact identity I5 computed at write time, files integrity-checked in run 2): u3/u5 of f4,5,7,9,10,12,14,15,16,17,52,53,54,56,58,59,61,63,64,65,66,67 (44 targets).
  **Verified in run 2**: u3.f68 and u5.f68, f69 ... f113 (30 targets). For run-1 targets W and a_j used by t05 were read from the integrity-checked core files, not from the old wvecs.pkl
  (not trusted after the reboot); `cert_summary.pkl` for those 44 is rebuilt from the store (frame_terms, counts), not from fresh recomputation.
- For run-1 targets the I5 identity was NOT recomputed in run 2; the evidence is run 1's log line plus xz CRC/sha256 and the independent smoke test on u3.f4.
  The independent full check of all 74 is left to the other checker.

## Wall times
t03 29 min (2+2 sectors in parallel); t04 run 1 ~2 h (45 targets); t04 run 2 75 min; t05 114 min. 4 threads, machine shared (load 15-140).

## Notes for a checker
- Row files are read from existing folders (W5C1, A3h4_frame_ds, A3h_weight4_target, A3h5_exact_certificate), see FORMAT.md; they are not copied into the store.
- sha256 values in MANIFEST.txt were computed after t05 finished; core files carry the t section as a second xz stream.
- Left on disk: `ysol_*.pkl` (3 GB, kept as instructed). Deleted: `ainv_*`, `wvecs.pkl`. Original W5G2 folder and paper/ untouched.
- Scripts changed relative to W5G2: g2common (HERE/SRC), t03 (reads trim/sel from W5G2), t04 (persist + RESUME), t05 (persist t_k); certio.py new.
