# W5G3R: exact re-verification of the 74 stored certificates (W5G3_cert_store)

"Refereed" here means an independent AI re-check with its own code. It is not peer review.

## Verdict

**CONFIRMED EXACTLY.** For all 74 targets v = u.f_i (u in {3,5}, 37 rows), my own reader and verifier,
written from `FORMAT.md` alone, check both identities in exact integer arithmetic over Z, after clearing all
denominators of each target. That is equivalent to checking them over Q.
- (I5) v = Σ a_j F_j + Σ c·act(r,g) + e ш W holds on all 46656 weight-5 words: **74/74, 0 mismatching words**.
- (I4) W = Σ t_k G_k holds on all 5184 weight-4 words: **74/74, 0 mismatching words**.

No modular arithmetic is used anywhere in this result. Both checks read the same W from the same core file.
This closes the weak point W5G2R recorded at L6, that the W of t04 and the W of t05 could not be compared on
stored data. Every row these identities reference was compared with an independent regeneration and found equal,
exactly and in full (no sampling). The standard rows and the frame rows were regenerated with the earlier
referee W4RX's code, and the lifts with my own code (§3). What remains taken on trust is listed in §5.

## 1. What I read and what I did not

- Read: `FORMAT.md`, `REPORT.md`, `MANIFEST.txt`; W5G2R's `REFEREE.md`, `rc.py`, `r02a_lifts.py`, `r03_w4.py`,
  `r04_frames.py`; W4RX's `result.md`, `w4rx_core.py`, `w4rx_rel.py`, `check_b.py`.
- Not read and not imported: `certio.py`, `t04_assemble.py`, `t05_w4.py`, `check_store_example.py`, and every
  other script of the authors. No module of the notes is imported. The pickles are read with a restricted
  unpickler that replaces the `qsqrt2.QS2` class with a data stub, so no project code runs when they are loaded.
- Python: project `.venv`. My own code is `common.py` (word enumeration in `itertools.product(range(-1,8))` order
  with first letter ≠ 0 and last letter ≠ −1, `act`, recursive shuffle product, primitive normalisation, xz/JSON
  reader, manifest check) and `verify.py` (the verifier). `regen_*.py` hold the regenerations and
  `controls.py` holds the controls.

## 2. The exact check (`verify.py`, log `verify.log`, one RESULT line per target)

For each target the script does the following.
1. It recomputes the sha256 of both store files and compares them with `MANIFEST.txt`. Result: 148/148 equal.
2. It parses every "p/q" as a pair of Python ints with q > 0. It checks the section structure, rejects duplicate
   terms, and checks the term counts against the meta line.
3. For (I5) it sets D = lcm of every denominator among a_j, c and W (7259 to 7264 digits). It accumulates
   D·(right side) on the 46656 words as Python integers:
   - literal rows: `rels_<family>.pkl[idx]`, mapped through the permutation of word indices for g = 7;
   - framew5 rows;
   - frame234 rows: their stored label is checked against the certificate's (n, rowidx, word);
   - e ш W, computed with my own general shuffle of (2) and (6) with each word of W.
   It then compares the result word by word with D·act(framew5 row i, u). Any word outside the 46656 admissible
   words would raise an error. None did.
4. For (I4) it sets D4 = lcm of the denominators of W and t (7270 to 7275 digits). It forms
   Σ (D4 t_k) G_k with G_k = T0[S[i]] (coefficients (0,c), q = 0 asserted) or A3h4 row ш word (my shuffle,
   not normalised, as FORMAT.md says). It compares the result with D4·W on all 5184 words.

| item | method | result |
|---|---|---|
| (I5), 74 targets | exact over Z (equivalent to over Q), all 46656 words | **0 mismatches, 74/74** |
| (I4), 74 targets | exact over Z, all 5184 words | **0 mismatches, 74/74** |
| sha256 of the 148 files | compared with MANIFEST.txt | 148/148 equal |
| g values used | read from the lit files | g ∈ {1, 7} only; std rows only with g = 1; lift1, lift2old, lift4old with g = 1 and g = 7; no u = 3 or u = 5 image of any row appears |
| terms per target | read from the core and lit files | framew5 59–64, frame234 13; std 38365–38367, lift1 8268, lift2old 193, lift4old 4721; W support 5180 of 5184 words (all nonzero); A3h5_S 5040, A3h4_frame 16–28 |
| lift words in frame234 / A3h4_frame terms (`check_twords.log`) | checked directly | all 29 distinct words admissible and of the right length |

Controls (`controls.log`; in memory, nothing on disk changed; certificate of u3.f4):

| control | I5 mismatching words | I4 mismatching words |
|---|---|---|
| unmodified | 0 | 0 |
| one std coefficient c (#1000) + 1 | 12 | 0 |
| one frame coefficient a_j (#3) + 1 | 2409 | 0 |
| one W coefficient (#17) + 1 | 10 | 1 |
| target replaced by its conjugate act(v,7) (= u5.f4) | 7236 | n/a |
| target replaced by another target, u3.f5 | 4496 | n/a |
| target replaced by the untwisted row f4 | 5476 | n/a |

## 3. Are the referenced rows the rows they claim to be? (all exact, all rows, no sampling)

| family (as used by the certificates) | how checked | result | status |
|---|---|---|---|
| std, weight 5 (94762 rows; regularised double shuffle + distribution d = 2, 4, 8) | regenerated with W4RX's generator `w4rx_rel.standard_rows(5, {})`, which imports only `w4rx_core` and no module of the notes (`regen_std.log`, 327 s) | identical ordered list, 94762/94762 | **regenerated independently and equal** (with an earlier referee's code, W4RX, read and understood; not my own) |
| lift1, lift2old, lift4old, weight 5 (20736 / 31104 / 35294) | own shuffle, primitive(b ш w) for every base b and admissible w (`regen_lifts.log`) | all equal label by label; the label sets are identical | **regenerated independently and equal** (own code) |
| lift bases B1 (4), B2 (54), B3 (539) | must be among W4RX's regenerated weight-1/2/3 rows (cumulative, as in W4RX check_b) | all contained | regenerated (W4RX code) and equal |
| lift base B4 = A3h5_S = T0[S] (5042 rows: std 4810, lift1 223, lift2old 7, lift3old 2); these are also the weight-4 generators G_k of (I4) | std rows from W4RX's `standard_rows(4, {})`; lifts with own shuffle | 5042/5042 literally equal | **regenerated independently and equal** |
| frame234 (all 32520 rows, of which 15 occur with nonzero a_j, on bases (2,3), (3,6)) | own shuffle, primitive(A3h4 row ш word) | 32520/32520 equal | lift step: own and equal; bases: next row |
| A3h4 frame rows, weights 2–4 (all 44 + 92 + 184 rows, including the 13 used in G_k and the 2 frame234 bases) | rebuilt from the stored provenance (frame and stuffle factors) with W4RX's conversion code: function-level double shuffle, Möbius/Cayley conversion, Q(√2); its convergence assertions are kept (`regen_frames.py`) | 320/320 equal; all convergence asserts pass; max stuffle vertex 0.8284 | **regenerated independently and equal** (W4RX code) |
| framew5 (all 288 rows: the 66 with nonzero a_j and the 37 target rows f_i) | same, W4RX conversion code, applied at weight 5 for the first time | 288/288 equal; all convergence asserts pass; max stuffle vertex 0.8284 | **regenerated independently and equal** (W4RX code) |
| A3h4_frame weight-4 lifts in G_k | computed by my verifier from the A3h4 rows (own shuffle) | used directly | own |

Comparison with the earlier referees: W5G2R regenerated the lifts (r02a). It checked the frame rows with the
authors' own conversion code (`frame_ds.convert`) and its own domain test (r04). It did not regenerate the std
rows; it only evaluated them numerically. My std and frame regenerations use W4RX's independent code instead,
so they are a stronger check. My used-row counts (66 framew5, 15 frame234 on bases (2,3), (3,6), the 13 A3h4
bases) agree with W5G2R's (`used_rows.log`).

## 4. What the check establishes and what it does not

It establishes three things:
1. The stored numbers satisfy (I5) and (I4) exactly.
2. Every row the identities reference equals a row produced by code independent of the authors'. The std rows and
   the frame conversion come from W4RX's code; the lifts, the shuffles and the verifier are my own.
3. Each frame row's stored provenance passes W4RX's convergence assertions.

Together with the standard arguments (W5G2R L4–L6, which I re-read and accept), this gives the following. Lifts by
admissible words are relations. Complex conjugation maps an integer relation on admissible words to a relation.
e ш W is a relation because W lies in the span of the G_k. Hence each u.f_i is a relation, given the
hypotheses in §5.

It does not establish:
- that the generated rows are true relations among the numbers H(w); that is a matter of theorems and of the
  correctness of the generating code, see §5;
- anything about the downstream use: W5C1's vector D, W5G1's split D = J + e ш Z + X, or W5C2's multipliers.

## 5. Taken on trust (not re-examined by me)

1. The cited theorems: regularised (extended) double shuffle and distribution relations at level 8 (Racinet,
   Zhao, Ihara–Kaneko–Zagier type), the shuffle product of convergent iterated integrals, and the stuffle product
   of absolutely convergent multiple polylogarithm series.
2. That the double shuffle/distribution specification shared by the authors' `relations.py` and W4RX's
   `w4rx_rel.py` is a correct instance of those theorems: the regularisation convention, the H ↔ G dictionary,
   and the signs. Two independent implementations agree on all 94762 + 8262 rows, but both follow one
   specification. W4RX states that it shares the enumeration order and normalisation with the notes.
3. The correctness of W4RX's frame conversion code (`w4rx_core`: Cayley map, path conversion, homotopy check) and
   of its convergence criterion. Two independent implementations (the authors' `frame_ds` and W4RX) agree on all
   608 frame rows. This replaces hypothesis (b), "correctness of the frame conversion code", by the agreement of
   two implementations; it does not prove either one correct.
4. The provenance data themselves (frames and factors) and the lower base selections, read as data. Each is
   checked only through the regenerated rows.
5. Numerical vanishing of the rows and targets: W5G2R r01, not rerun. It is evidence, not needed for the logic.
6. The standard arguments for conjugates and for e ш (relation) (W5G2R L5, L6). I re-read them and accept them,
   but I did not formalise them.
7. Python's integer arithmetic and its lzma, json and hashlib modules.
8. Everything downstream of "each u.f_i is a relation": W5C1's D, W5G1's exact split (independently checked only
   modulo two primes, by W5G1R), and W5C2's multipliers.

## 6. Timings

- One target: about 40 s (read and parse 6 s, literal rows 30 s, rest 4 s), plus about 3 s for the sha256 values.
  Peak memory is 5.2 GB per process.
- Full run: 3 worker processes, 15:41:35–15:58:28 (17 min wall clock, 2982 s CPU in total). I used at most
  4 cores including the regenerations.
- Regenerations: std rows at weight 5 327 s; lifts and T0[S] 28 s; frame rows 5 s; used rows 110 s; controls 4.5 min.
- Wall time of the whole referee task: about 25 min, 15:36–16:01.
- The folder is 2.4 MB. No temporary files were left, and nothing outside this folder was written.

## 7. Sentence for the paper

To replace "it did not repeat the exact computation (refereed)":

> A second independent re-check repeated the exact computation with its own code: it verified the stored
> certificates of all $74$ images in exact rational arithmetic on all $46656$ weight-five and all $5184$ weight-four
> words, after regenerating every relation they use, by code independent of ours, and finding it identical (refereed).

(This sentence does not remove the standing hypotheses. The cited theorems and the correctness of the
relation-generating and conversion code, now implemented twice, remain assumptions.)

## Files

`common.py`, `verify.py` + `run_verify.sh` (resumable; skips targets that already have a RESULT line) → `verify.log`;
`controls.py` → `controls.log`; `regen_std.py` → `regen_std.log`, `regen_std.pkl` (2 MB, W4RX rows of weights 1–4);
`regen_lifts.py` → `regen_lifts.log`; `regen_frames.py` → `regen_frames.log`, `regen_frames.pkl`;
`used_rows.py` → `used_rows.log`, `used_rows.json`; `check_twords.py` → `check_twords.log`.
