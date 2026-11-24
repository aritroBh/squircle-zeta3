# W4RX: independent re-implementation of the certificate of thm:W4R (W_4 = R)

**VERDICT: INDEPENDENTLY CONFIRMED.** All 5068 certificate rows and both target vectors were rebuilt by
code that imports no module of the notes; they agree exactly with the notes' rows and targets, and the
certificate identity holds with our own reading of `D_certificate.json` and our own exact arithmetic
(residual 0 on all 5184 words, both parts, no pi^4 row). Every converted quantity agrees with a direct
numerical value to about 147 digits, and every row vanishes numerically to 6e-151.

Run everything: `./run_all.sh` (about 50 s). Python: `../../../.venv/bin/python`.

## Coverage

| item | status | result (log) |
|---|---|---|
| GPL expression of W_4 (148 terms) | re-derived from the definition of W_4 (l1, l2, shuffles) | equal to `A3_W4_proof/gpl_expression.tsv` (check_a.log) |
| Cayley map, conversion lemma, vectors W, R, D (rational and sqrt2 parts) | re-implemented (`w4rx_core.py`, Goncharov G notation, H only for comparison) | equal word by word to `A3c.../Q/vectors.json` and to the verifier's rebuilt target (check_a.log) |
| conversion of each Li_n(z) | re-implemented; all admissible frames converted (12 to 68 per atom) | all frames agree numerically with mpmath to 7e-151; the frame of the notes (first in a fixed order) gives exactly the notes' vector |
| numerical values of all 5831 admissible level-8 words of weight <= 4 | re-implemented (Hoelder convolution at 1/2, python-flint acb, 600/660 bits, N = 900/1100) | the two runs agree to < 1e-150 (check_a.log) |
| W_4 numerically | tanh-sinh quadrature at 150 digits | = converted W vector to 4.9e-150 |
| each term of R (11 polylog atoms, L, lambda, pi^2, zeta(3), 29 products) | direct mpmath values | = converted vectors to <= 2e-149; R to 1e-147; W_4 - R(direct) = 2e-149 |
| standard rows weight 1, 2, 3 and all 8262 of weight 4 (DS + distribution) | re-implemented (`w4rx_rel.py`) | ordered lists identical to the notes' (check_b.log) |
| lower bases (4 rows weight 1, 54 rows weight 2) | read as data | every row is one of our regenerated relations |
| 230 lifts (223 lift1, 7 lift2old) | re-implemented | equal |
| 28 frame rows | re-derived from (frame, factors) read as data via a restricted unpickler; own stuffle/shuffle/Moebius/Cayley code; convergence checks (all z real, 0<|z|<1, every stuffle vertex real with |b|<1 and 1/b a frame point, no pole with t in (0,1)) | equal to the stored rows and to the certificate rows |
| certificate identity | own reading of D_certificate.json, own exact Fraction arithmetic, our rows and our target | residual 0 on all 5184 words in both parts; heights <= 381 digits; one zero multiplier (check_cd.log) |
| (d) numerics of every row | all 5068 certificate rows; also all 8262 standard rows | max |val| 5.9e-151 (rows), 1.4e-148 (all std); negative control (one coefficient +1) gives 0.02..1.08 |
| multipliers y_i | data, not recomputed | not needed for the proof |
| original verifier `verify_D.py` | RERUN (log redirected, no file written) | passes (rerun_original_verifier.log) |

Rows used, by the theorem they need (classify_rows.log):
- double shuffle, both factors convergent: 4032 (u of weight 1) + 100 (u of weight 2): elementary (proved in the section)
- double shuffle with divergent u: 576 (u=(0)) + 42 (u in {(0,0),(0,1),(0,2),(0,3),(0,4)}): **cited** regularised double shuffle (IKZ / Racinet / Zhao)
- distribution d=2: 35 convergent + 25 regularised: proved in the section (the regularised case is a new short argument)
- lift1 223 (basis rows: 2 weight-one log identities, 2 weight-one distribution), lift2old 7 (one weight-2 distribution row): proved
- frame 28 (10+9+9, from 1+3+9 distinct relations, frames (inf,1,-a) and (inf,a,b)): proved (function-level double shuffle, |z|<1)

## What is re-implemented, sampled, rerun

- **Re-implemented fully:** the conversion (W, all R atoms, constants), the numerical evaluator, all relation
  families that occur (weight-one, DS, distribution incl. regularisation, lifts, frame DS), the certificate check.
- **Shared with the notes (not independent):** the mathematics; the enumeration order and the row normalisation
  of the standard rows (a specification needed to read the labels; also checked by content); the choice of frame
  for each Li_n(z) (other frames give other, numerically equal vectors; with our own choice the vectors would
  differ by relations we did not try to certify); the data read: labels and multipliers, the frame provenance,
  the two lower bases.
- **Only rerun:** the notes' verifier and the notes' row builder (`extract_original_rows.py`, comparison only).
- **Sampled:** nothing; every row was rebuilt and evaluated.

## Files

- `w4rx_core.py` (Q(sqrt2), Cayley map, conversion lemma, shuffle, atoms, H/G dictionary), `w4rx_rel.py` (relations)
- `evalG.py` (numerics), `check_a.py` (part a), `check_b.py` (part b), `check_cd.py` (parts c, d),
  `classify_rows.py`, `frame_table.py`, `reg_identity_check.py`
- RERUN scripts (import the notes' code; comparison only): `extract_original_rows.py`, `rerun_original_verifier.py`
- `logs/`, `data/` (values, comparison rows, pickles; all < 5 MB)

## Notes for the referee

1. Own bug found and fixed during the work: the weight-one row for j=4 was built from a dict literal with a
   repeated key; the comparison with the notes' ordered weight-1 list caught it. Not in the notes' code.
2. The conversion needs no regularisation for any term used (outermost form regular at the end point, innermost at
   the start; the code asserts it). The PLANNER_AUDIT worry about regularised radial pieces does not arise.
3. The regularised double shuffle theorem (618 rows) is the only cited input; locators in Racinet 2002 /
   Zhao 2010 / IKZ 2006 not checked (no access).
4. Supplement section: `paper/ieee/proofs/sec_g1_W4R.tex`, wrapper `wrap_g1r.tex`, macros `macros_g1r.tex`,
   PDF `wrap_g1r.pdf` (13 pages incl. title and references; clean build).
