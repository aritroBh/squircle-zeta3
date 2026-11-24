# W5G2: exact certificates for the 74 conjugated frame relations

## (1) Verdict: **CLOSED**

All 74 vectors in (H) have exact rational certificates. Each was checked exactly, with no prime and no numerics in the check.
Each certificate uses only literally generated proved relations, plus iπ ш (weight-4 vector).
Every weight-4 vector that appears lies exactly in the span of proved weight-4 relations, so its iπ ш part is a relation too.
No Galois or transport argument enters any certificate.
The only group action used is u = 7, which is complex conjugation (W5G1 claim B, confirmed by the referee).

## (2) What is now proved

Notation:
- f_i is one of the 37 open frame-w5 rows, indices into `W5C1/validated_w5.pkl`:
  4, 5, 7, 9, 10, 12, 14, 15, 16, 17, 52, 53, 54, 56, 58, 59, 61, 63–75, 87, 91, 95, 103, 104, 107, 113.
- u.f_i, for u = 3, 5, is the image under the letter map j ↦ uj mod 8.
- e = H(2) − H(6). Its value is iπ/2.

**Theorem (G2).** For each of the 74 pairs (i, u) there is an exact identity on all 46656 admissible weight-5 words:

  u.f_i = Σ_j a_j F_j + Σ c_r r + Σ c'_r (7.r) + e ш W,   with   W = Σ_k t_k G_k  exactly on all 5184 weight-4 words.

The terms are:
- **F_j**: frame relations taken from a basis of 90:
  - 67 weight-5 frame rows from `validated_w5.pkl`;
  - 23 weight-5 lifts (frame row ш word) of A3h4 frame rows of weights 2–4, from `rels_frame234.pkl`.
  - In total the certificates use 66 distinct frame-w5 rows and 13 A3h4 frame rows (weights 2–4).
  - Each one was regenerated from its stored provenance by `frame_ds.candidates/identity/convert`, which enforces the domain conditions, and equals the stored row (t06).
- **r**: literal rows of the std family (weight-5 regularised double shuffle and distribution) and of the lift families `lift1`, `lift2old`, `lift4old` (proved base row ш word).
  - Every u = 3, 5 image of a std row was matched to a literal std row.
  - Every u = 3, 5 image of a lift row was rewritten, by an exact decomposition of its base row's image, as literal lift rows plus an e ш (·) part.
  - Those e-parts go into W.
- **7.r**: complex conjugates of literal lift rows. These are relations by W5G1 claim B.
- **G_k**: proved weight-4 relations:
  - A3h5's 5042 base rows T0[S] (Theorem W4R's inputs);
  - weight-4 lifts of A3h4 frame rows (2,3), (3,3), (3,5), (3,6), (4,4), (4,6), (4,8), (4,9), (4,11), (4,12), (4,40), (4,41), (4,43). These are the 13 A3h4 rows above, each regenerated in its domain.
  - Then e ш W = Σ t_k (e ш G_k). Each e ш G_k is a relation: all words are admissible, so its value is (iπ/2)·0.

**Corollary.** Statement (H) holds: every u.f_i is a relation.
W5G1's residual X is a rational combination of exactly these 74 vectors (W5G1 s05, re-derived by the referee). So X·V = 0.
Together with W5G1's exact split D = J + e ш Z + X and its certificate for Z, this gives D·V = 0, i.e. **W₅ = N₅ + R₅**.

The hypotheses are the same standing ones as Theorem W4R:
- (a) The std, distribution and regularisation families are theorems (Ihara–Kaneko–Zagier double shuffle; Zhao's distribution relations). Lifts of proved relations are relations.
- (b) The frame conversion code (`frame_ds.convert`, `cmzv.path`, Q(√2) arithmetic) is correct. This is W4R's caveat (iii). It now covers 66 weight-5 frame rows and 13 weight-2–4 frame rows, all regenerated in domain.
- (c) W5C1's exact vector D (conversion of W₅, N₅ and R₅, including the path-reversal duality for Li_{3,2}(1/√2, 1)) is correct. This is outside the Galois question and unchanged.
- (d) W5G1's exact identity D = J + e ш Z + X (s05, 0 mismatches). I did not re-run it here; the referee checked it modulo two primes.

There is no hypothesis about Galois images any more.

## (3) The 74 certificates

Column meanings:
- *frame terms*: number of nonzero a_j.
- *literal rows*: distinct proved rows r with nonzero coefficient, split as std / lift1 / lift2old / lift4old. The u = 7 conjugates are counted in the next column.
- *max coef digits*: digits of numerator plus denominator of the largest coefficient on literal rows.

The a_j themselves have at most 27 digits.
The last column gives the mismatches in the exact weight-5 identity and in the exact weight-4 identity W = Σ t_k G_k.

| # | vector | frame terms | literal rows (std/lift1/lift2old/lift4old) | of which u=7 conj. | max coef digits | W support | W: proved w4 rows (frame lifts) | W max digits | exact checks |
|---|---|---|---|---|---|---|---|---|---|
| 1 | u3.f4 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13835 | 5180 | 5056 (16) | 14558 | 0 / 0 mismatches |
| 2 | u5.f4 | 72 | 51548 (38366/8268/193/4721) | 1766 | 13841 | 5180 | 5056 (16) | 14562 | 0 / 0 mismatches |
| 3 | u3.f5 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5056 (16) | 14556 | 0 / 0 mismatches |
| 4 | u5.f5 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5056 (16) | 14558 | 0 / 0 mismatches |
| 5 | u3.f7 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14566 | 0 / 0 mismatches |
| 6 | u5.f7 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14566 | 0 / 0 mismatches |
| 7 | u3.f9 | 72 | 51548 (38366/8268/193/4721) | 1766 | 13842 | 5180 | 5056 (16) | 14561 | 0 / 0 mismatches |
| 8 | u5.f9 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5056 (16) | 14564 | 0 / 0 mismatches |
| 9 | u3.f10 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13843 | 5180 | 5056 (16) | 14560 | 0 / 0 mismatches |
| 10 | u5.f10 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13843 | 5180 | 5056 (16) | 14562 | 0 / 0 mismatches |
| 11 | u3.f12 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13843 | 5180 | 5056 (16) | 14564 | 0 / 0 mismatches |
| 12 | u5.f12 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13841 | 5180 | 5056 (16) | 14563 | 0 / 0 mismatches |
| 13 | u3.f14 | 72 | 51548 (38366/8268/193/4721) | 1766 | 13839 | 5180 | 5056 (16) | 14562 | 0 / 0 mismatches |
| 14 | u5.f14 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13840 | 5180 | 5056 (16) | 14564 | 0 / 0 mismatches |
| 15 | u3.f15 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13841 | 5180 | 5056 (16) | 14562 | 0 / 0 mismatches |
| 16 | u5.f15 | 72 | 51548 (38366/8268/193/4721) | 1766 | 13843 | 5180 | 5056 (16) | 14561 | 0 / 0 mismatches |
| 17 | u3.f16 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13840 | 5180 | 5056 (16) | 14562 | 0 / 0 mismatches |
| 18 | u5.f16 | 72 | 51548 (38366/8268/193/4721) | 1766 | 13840 | 5180 | 5056 (16) | 14562 | 0 / 0 mismatches |
| 19 | u3.f17 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13843 | 5180 | 5056 (16) | 14564 | 0 / 0 mismatches |
| 20 | u5.f17 | 72 | 51549 (38367/8268/193/4721) | 1766 | 13839 | 5180 | 5056 (16) | 14564 | 0 / 0 mismatches |
| 21 | u3.f52 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14560 | 0 / 0 mismatches |
| 22 | u5.f52 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14562 | 0 / 0 mismatches |
| 23 | u3.f53 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14560 | 0 / 0 mismatches |
| 24 | u5.f53 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14563 | 0 / 0 mismatches |
| 25 | u3.f54 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14565 | 0 / 0 mismatches |
| 26 | u5.f54 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14560 | 0 / 0 mismatches |
| 27 | u3.f56 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14561 | 0 / 0 mismatches |
| 28 | u5.f56 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14563 | 0 / 0 mismatches |
| 29 | u3.f58 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13843 | 5180 | 5065 (25) | 14563 | 0 / 0 mismatches |
| 30 | u5.f58 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13843 | 5180 | 5065 (25) | 14563 | 0 / 0 mismatches |
| 31 | u3.f59 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14563 | 0 / 0 mismatches |
| 32 | u5.f59 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13843 | 5180 | 5065 (25) | 14565 | 0 / 0 mismatches |
| 33 | u3.f61 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14564 | 0 / 0 mismatches |
| 34 | u5.f61 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14562 | 0 / 0 mismatches |
| 35 | u3.f63 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14564 | 0 / 0 mismatches |
| 36 | u5.f63 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14562 | 0 / 0 mismatches |
| 37 | u3.f64 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14560 | 0 / 0 mismatches |
| 38 | u5.f64 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14564 | 0 / 0 mismatches |
| 39 | u3.f65 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13838 | 5180 | 5065 (25) | 14559 | 0 / 0 mismatches |
| 40 | u5.f65 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13838 | 5180 | 5065 (25) | 14561 | 0 / 0 mismatches |
| 41 | u3.f66 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13838 | 5180 | 5065 (25) | 14561 | 0 / 0 mismatches |
| 42 | u5.f66 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13838 | 5180 | 5065 (25) | 14563 | 0 / 0 mismatches |
| 43 | u3.f67 | 75 | 51547 (38365/8268/193/4721) | 1766 | 13843 | 5180 | 5067 (27) | 14567 | 0 / 0 mismatches |
| 44 | u5.f67 | 75 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5068 (28) | 14567 | 0 / 0 mismatches |
| 45 | u3.f68 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13837 | 5180 | 5065 (25) | 14560 | 0 / 0 mismatches |
| 46 | u5.f68 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13841 | 5180 | 5065 (25) | 14562 | 0 / 0 mismatches |
| 47 | u3.f69 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14567 | 0 / 0 mismatches |
| 48 | u5.f69 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14567 | 0 / 0 mismatches |
| 49 | u3.f70 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13843 | 5180 | 5056 (16) | 14564 | 0 / 0 mismatches |
| 50 | u5.f70 | 72 | 51549 (38367/8268/193/4721) | 1766 | 13843 | 5180 | 5056 (16) | 14564 | 0 / 0 mismatches |
| 51 | u3.f71 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14566 | 0 / 0 mismatches |
| 52 | u5.f71 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13837 | 5180 | 5065 (25) | 14560 | 0 / 0 mismatches |
| 53 | u3.f72 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14566 | 0 / 0 mismatches |
| 54 | u5.f72 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13841 | 5180 | 5065 (25) | 14566 | 0 / 0 mismatches |
| 55 | u3.f73 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14566 | 0 / 0 mismatches |
| 56 | u5.f73 | 72 | 51548 (38366/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14566 | 0 / 0 mismatches |
| 57 | u3.f74 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13836 | 5180 | 5056 (16) | 14560 | 0 / 0 mismatches |
| 58 | u5.f74 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13836 | 5180 | 5056 (16) | 14560 | 0 / 0 mismatches |
| 59 | u3.f75 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13838 | 5180 | 5065 (25) | 14558 | 0 / 0 mismatches |
| 60 | u5.f75 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13838 | 5180 | 5065 (25) | 14561 | 0 / 0 mismatches |
| 61 | u3.f87 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14567 | 0 / 0 mismatches |
| 62 | u5.f87 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14567 | 0 / 0 mismatches |
| 63 | u3.f91 | 72 | 51547 (38365/8268/193/4721) | 1766 | 13838 | 5180 | 5056 (16) | 14564 | 0 / 0 mismatches |
| 64 | u5.f91 | 72 | 51549 (38367/8268/193/4721) | 1766 | 13841 | 5180 | 5056 (16) | 14564 | 0 / 0 mismatches |
| 65 | u3.f95 | 73 | 51547 (38365/8268/193/4721) | 1766 | 13837 | 5180 | 5065 (25) | 14565 | 0 / 0 mismatches |
| 66 | u5.f95 | 73 | 51548 (38366/8268/193/4721) | 1766 | 13829 | 5180 | 5065 (25) | 14562 | 0 / 0 mismatches |
| 67 | u3.f103 | 77 | 51547 (38365/8268/193/4721) | 1766 | 13841 | 5180 | 5067 (27) | 14564 | 0 / 0 mismatches |
| 68 | u5.f103 | 77 | 51547 (38365/8268/193/4721) | 1766 | 13841 | 5180 | 5068 (28) | 14564 | 0 / 0 mismatches |
| 69 | u3.f104 | 74 | 51547 (38365/8268/193/4721) | 1766 | 13840 | 5180 | 5065 (25) | 14567 | 0 / 0 mismatches |
| 70 | u5.f104 | 74 | 51548 (38366/8268/193/4721) | 1766 | 13838 | 5180 | 5065 (25) | 14563 | 0 / 0 mismatches |
| 71 | u3.f107 | 72 | 51548 (38366/8268/193/4721) | 1766 | 13842 | 5180 | 5065 (25) | 14564 | 0 / 0 mismatches |
| 72 | u5.f107 | 72 | 51548 (38366/8268/193/4721) | 1766 | 13839 | 5180 | 5065 (25) | 14562 | 0 / 0 mismatches |
| 73 | u3.f113 | 81 | 51547 (38365/8268/193/4721) | 1766 | 13844 | 5180 | 5059 (19) | 14565 | 0 / 0 mismatches |
| 74 | u5.f113 | 81 | 51547 (38365/8268/193/4721) | 1766 | 13844 | 5180 | 5059 (19) | 14565 | 0 / 0 mismatches |

All 74 pass: weight-5 identity 0 mismatches (t04.log), weight-4 identity 0 mismatches (t05.log).
Certificates for u = 3 and u = 5 of the same row are nearly identical, as expected: Fold_s(u.f) = χ_s(u) Fold_s(f).

## (4) How it was done and how to rerun

Python: `../../../.venv/bin/python`. Run in this folder, in this order.

1. `t01_trim.py pp|pm|mp|mm` (4 in parallel, about 52 min).
   - Folds the G-stable generating set Et = std + lift1 + I (I = e ш words(4)) + lift2old/3old/4old into each sector.
   - Runs sparse echelon mod p with pivot rows (W5C2's LeanEchelon), then builds the null basis of the quotient.
   - Ranks: 12120 / 11065 / 11057 / 12006, nullity 408 in total, the same as W5G1 s08.
   - The null basis was checked against every folded row.
   - Output: `trim_*.pkl`, `quot_*.npz`.
2. `t02_select.py` (2 min).
   - Works mod p. Picks a basis of 90 for the frame candidates in the quotient: 288 frame-w5 rows, 32520 frame-234 lifts and their u = 7 images. The basis has 67 frame-w5 rows and 23 frame-234 lifts.
   - All 74 targets lie in its span. This only guides the choice of basis.
   - Output: `sel.pkl`.
3. `t03_solve.py <s>` (two at a time, about 12 min each).
   - Exact Dixon lifting (`dixon2.py`: flint inverse mod p, float-BLAS lifting, common-denominator reconstruction, exact integer check M·n = d·B).
   - Right-hand sides: the 37 rows f_i and the frame basis, 90 in all. Fold_s(u.f) = χ_s(u) Fold_s(f), so these are enough.
   - Gives exact y and exact normal forms on the free columns. Steps 496–512; heights about 1750 digits.
   - Output: `ysol_*.pkl`, about 750 MB each (deleted).
4. `t04_assemble.py` (about 100 min).
   - Exact quotient solve: NF(u.f) = Σ a_j NF(F_j) over Q. The rank is 90, the a_j are found and checked exactly, and they have at most 27 digits.
   - Unfolds x = u.f − Σ a_j F_j as Σ_s P_s x and rewrites every g.r as literal proved rows.
   - The base-row decompositions are exact: weight 1 (6 pairs; the e-part is nonzero only for the 2 branch rows), weight 2 (4), weight 4 (6).
   - **Verifies each identity exactly on all 46656 words** (integers after one common denominator).
   - Output: `cert_summary.pkl`, `wvecs.pkl` (1.8 GB, deleted).
5. `t05_w4.py` (about 100 min).
   - Uses 8686 weight-4 generators, rank 5070.
   - Dixon on the 5070 × 5070 square system with all 74 W, then **an exact check of W = Σ t_k G_k on all 5184 words**.
   - Output: `w4cert_summary.pkl`.
6. `t06_frames.py` (30 s).
   - Regenerates every frame row used in t04/t05 from its provenance, with W4R verify_D's domain assertion (`fs in candidates(n, bs, arity)`), and compares it with the stored row.
   - 66 + 13 rows, all equal; max |z| = 0.83.

Helpers:
- `g2common.py` and `g2lin.py`: `bases` and `exact_combination` copied from W5G1.
- `w5c2_fold.py` and `w5c2_lib.py`: copied from W5C2 unchanged.

Every step logs to `*.log` and `*.out`. The large intermediates were deleted (ainv_*.npy, ysol_*.pkl, wvecs.pkl); steps 3–5 regenerate them. The total wall time is about 5 h on this machine (14 cores, 48 GB).

Notes on rigour:
- The mod-p steps (t01, t02) only choose rows, columns and a basis. Every claim rests on the exact checks in t03 (M·n = d·B over Z), t04 (the full weight-5 identity over Q) and t05 (the full weight-4 identity over Q).
- The sector fold and unfold is used only to find the certificate. The verified identity is a plain sum of literal rows on unfolded words.

## (5) What a paper sentence may now claim

- The images of the 37 frame relations under √2 ↦ −√2 are relations. We prove this with exact computer certificates.
- Each image equals a rational combination of proved weight-five relations, plus iπ times a shuffle with proved weight-four relations.
- No Galois argument is used. Only complex conjugation is used, and it maps relations to relations.
- So the closed form of W₅ with three extra arguments is proved, computer-assisted.
- It rests on the same assumptions as Theorem W4R: the cited double-shuffle and distribution theorems, and the correctness of our frame conversion code.
- It also rests on the exact conversion of W₅, N₅ and R₅ to words.

Suggested replacement for the "conditional" sentence:
"The sector decomposition also uses the images of the relations under √2 ↦ −√2. These images need not be relations, so we certified each one used directly. Each is an exact rational combination of proved weight-five relations and of iπ times shuffles of proved weight-four relations, checked word by word. Hence the closed form of W₅ is computer-assisted and holds under the same standing assumptions as Theorem W4R."

The open problem "prove that the images of the 37 weight-five frame relations … are relations" can be removed.
