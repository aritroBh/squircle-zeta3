# W5C14: can the fourteen-argument form of W_5 be proved?

## Verdict: PROVED (computer-assisted), with the same standing hypotheses as W_5 = N_5 + R_5

- R_14 − R_5 = 0 holds as an exact identity of level-8 periods. The certificate has the same structure as the one for R_5 (W5C2 + W5G1 + W5G2).
- With the proved W_5 = N_5 + R_5, this gives **W_5 = N_5 + R_14**, the fourteen-argument form displayed in the paper (main.tex, ~l. 995).
- Nothing new is assumed. The hypotheses are W5G2R's (a)–(d): the cited double shuffle/distribution theorems and the row code, the frame conversion code, and the exact conversion of atoms to words. In addition there are three new atom conversions; each one is a change of variables in an iterated integral (below).

## The statement

- R_14 = `M4f_w5_extras/R5_best.json` (70 classical terms). This is the paper's form.
  - Its Li_5 arguments are the 7 of R, plus 2b, (2−√2)/4, (2+√2)/4, −b/2, (4−3√2)/8, b⁴ and 1−b⁴.
  - Every ≥4-digit number in the JSON appears in `paper/ieee/supplement_W5.tex` (71/71; the copy is `paper_supplement_W5.copy.tex`).
- R_5 = `M4f_w5_extras/R5_min3.json` (66 terms). This is the form certified in W5C1/W5C2/W5G2.
- R_14 − R_5 has 28 nonzero classical terms (listed in s04.log). Among them are Li_5 at 14 arguments, ζ(5) and L(χ_8,5). (Corrected 2026-11-23 from 13; s04.log lists 14; found by the referee W5C14R.)

## Step 1: data and numerics (numerical)

- s03 (360 digits): N_5+R_14−W_5 = −6.6e−387, N_5+R_5−W_5 = −1.1e−386, R_14−R_5 = 4.6e−387.
- W_5 is from `M4f_w5_extras/W5_5400.json`, and N_5 is from `M4e_w5_lchi/cobracket.json`.

## Step 2: which arguments are level-8 (exact computation)

- s01: there are 27 distinct cross-ratios of four of the 8 points. The check over all ordered 4-tuples was exhaustive.
- **Cross-ratios (in):** b, −b, −b², 1/√2, 1/2, 1−1/√2, 2−√2, −b/2, (2∓√2)/4, **2b** (e.g. (0,a,b,∞)), b², −1, −1/√2.
- **Not cross-ratios:** (4−3√2)/8, b⁴, 1−b⁴.
- Why the three extras were certifiable: they are cross-ratios, so the converter `cmzv.find_cross` handles them. 2b is also a cross-ratio.
- The other three failed only because the converter knew cross-ratios and nothing else. This is not a real obstruction:
  - Let z = b⁴. Then (4−3√2)/8 = z/(z−1), and 1−b⁴ = 1−z.
  - Pull back along t = (v/a)². This sends {0, ±1, ∞, b²} into the 8 points (v = ±a, and b² ↦ b). s02 finds this map and 15 others.
  - So each of the three is an iterated integral on P¹ minus {0, ±a, ∞}, between points of the configuration, along a real segment that meets no letter:
    - Li_5(b⁴) = I_{0→b}[(2 dlog v)^4 ∘ −dlog((v−a)(v+a))];
    - Li_5((4−3√2)/8) = I_{0→b}[(2 dlog v − dlog(v−a) − dlog(v+a))^4 ∘ (dlog(v−a)+dlog(v+a))], via t = s/(s−1), s = (v/a)²;
    - Li_5(1−b⁴) = I_{a→b}[(dlog(v−a)+dlog(v+a))^4 ∘ −2 dlog v], via t = 1−s. The segment [b,a] contains 1, but 1 is not a letter.
  - These identities are proved by the change of variables in iterated integrals. The Li_5(b⁴) case is the duplication relation, applied twice.
  - Numerically (s05, arb): the three conversions match mpmath to 3e−165. The whole vector Δ evaluates to 3e−161.

## Step 3: the certificate (exact computation, all exit 0)

1. **s04:** Δ = R_14 − R_5 as exact rational and √2 vectors on admissible weight-5 words, using the project's `cmzv` conversion.
   - Supports are 12181/12181. All words are admissible.
2. **s06–s07:** Δ is folded into the four Galois sectors with W5C2's fold code.
   - The linear system is the one already trimmed in W5C2: the same rows S and pivot columns C. These were read from `W5C2_exact_certificate/work`; nothing was re-trimmed.
   - The solve is a Dixon solve with both components at once. Heights are up to 2080 digits. The check A·Y = B mod a second prime passed.
   - Time: 13–20 min per sector.
3. **s08:** Σ y_j R_{S_j} = Φ_s(Δ) **exactly** (Fractions) on all sector columns. All four sectors have 0 mismatching columns.
   - Every used sector row was recomputed as (1/λ)·fold of its unfolded source row (W5C1 rels_*.pkl / validated_w5.pkl).
   - Frame-w5 rows used: 29 (pp), 31 (pm), 0 (mp), 0 (mm). They are all among W5C2's 38.
4. **s09:** W5G1's s05 split was applied to Δ, using W5G1's per-base-row `lift_decomp.pkl` (read-only).
   - Exact on all words, both components, 0 mismatches: Δ = (literal std/lift/frame rows) + (u=7 conjugates of literal rows) + e ш Z′ + X′.
   - u=3,5 images of std rows: all literal std rows. Images of frame row 80: literal.
5. **s10:**
   - X′ is a combination of u.f_i for 30 frame rows × u ∈ {3,5} = 60 pairs. **All 60 are among the 74 pairs certified exactly in W5G2.** 0 pairs are uncertified.
   - Z′ (weight 4, support 646) is an **exact** combination of 5047 proved weight-4 rows: A3h5's T0[S] plus lifts of the A3h4 frame rows (2,3), (3,3), (3,5).
   - So e ш Z′ is a relation, by the same argument as W5G2R L6.

Conclusion: every term is a proved relation (W5G2R links L3–L7), so Δ·V = 0. That is R_14 = R_5, hence W_5 = N_5 + R_14. **Proved, computer-assisted.**

## What it rests on

- Cited theorems: regularised double shuffle and distribution at level 8 (Racinet; IKZ-type; Zhao), shuffle of convergent iterated integrals, stuffle/shuffle of absolutely convergent series, and change of variables in iterated integrals.
- Code correctness, unchanged from W4R/W5G2R:
  - `relations.py` rows;
  - `frame_ds`/`cmzv` conversion, now also used by my three pullback forms (s04, NEW dict);
  - W5G1's lift_decomp and W5G2's 74 certificates. Their exactness is as stated in their refereed results; W5G2's W was not stored (W5G2R caveat).
- Not independently refereed yet.

## Files and rerun (cwd = this folder, PY = ../../../.venv/bin/python)

- `s01_crossratio.py` (cross-ratio table)
- `s02_embed.py` (Möbius embeddings for the pullback)
- `s03_numeric.py 360`
- `s04_convert.py` (→ delta.pkl, 40 s)
- `s05_validate.py` (arb check, 30 s)
- `s06_fold.py`
- `run_s07.sh` (→ s07_dixon.py per sector, ~1 h total, ~5 GB)
- `s08_verify.py` (→ used_rows_delta_{sec}.pkl, 1 min)
- `s09_split.py` (24 min, writes the 750 MB split_delta.pkl)
- `s10_Zcert.py` (3 min)
- Helpers: `common.py`, `qs2.py`.
- Logs: s0N*.log.
- Kept: delta.pkl, delta_folds.pkl, yexact_{pp,pm,mp,mm}.pkl (the sector multipliers), Zprime_certificate.pkl.
- Deleted after the run (regenerated by s08/s09): used_rows_delta_*.pkl and split_delta.pkl.
- Nothing outside this folder was written.
