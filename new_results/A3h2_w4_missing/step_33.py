WORK='<project>/research/08_theory/A3h2_w4_missing'
md='''# A3h2 — the 29 missing weight-4 level-8 relations and the D-target (IN PROGRESS)

Workdir: `research/08_theory/A3h2_w4_missing`. Progress log: `progress.log`.

## Headline so far
- **Galois trick (brief §"The Galois trick") is only HALF true.** Proved form-level engine rows
  (shuffle/distribution/branch) are invariant under all σ_k, k∈{1,3,5,7}, but the Rogers five-term
  row and R1–R5 vanish only at k=1 and k=7 (complex conjugation); at k=3,5 they fail by explicit
  lower-weight periods. Consequence: the 8-constraint Galois LLL can only find the FORM-LEVEL
  subset of the missing relations; the full search must use single-embedding (k=1) LLL at high
  precision. Details below.
- Evaluator (A3h Q/check_conversion engine, traps reset) validated: matches A3f4's certified
  `vals_hi.pickle` (750 digits) to 5.6e-151 on the 32 weight-3 free words; ζ(3) anchor to 1.5e-100+;
  depth-4 word (7,4,7,7) confirmed to 1.7e-38 by an independent nested `mp.quad` (~10 min).
- 110 free words of the T3 quotient evaluated at all 4 embeddings, N=2000 (certified vs N=1400 to
  ≥ 2.2e-264 worst case, extrapolated ~369 digits). Saved: `vals110.pickle`.
- LLL pipeline validated on the known weight-3 case: flint fmpz_mat.lll with C=10^170 recovers the
  5-dim R1–R5 lattice with a 1e17 gap to the noise floor.

## Step 1 — Galois trick validation (100-digit target; prec 700 bits, N=800)
Setup: `evsetup.py` resets ctx.prec/cap, recomputes roots, clears caches (the A3f3/A3f4 traps).
σ_k on words: letter j ↦ k·j mod 8, −1 fixed; on coefficients a+b√2: √2 ↦ ±√2 with sign
GAL = (k=1:+, 3:−, 5:−, 7:+).

| test | k=1 | k=3 | k=5 | k=7 |
|---|---|---|---|---|
| Rogers five-term (36 words) | 1.4e-120 | **11.08** | **11.08** | 1.4e-120 |
| R1 | 8.8e-119 | **22.8** | **22.8** | 8.8e-119 |
| R2 | 3.2e-118 | **185** | **185** | 3.2e-118 |
| R3 | 3.3e-118 | **178** | **178** | 3.3e-118 |
| R4 | 5.4e-118 | **42.2** | **42.2** | 5.4e-118 |
| R5 | 8.8e-118 | **250** | **250** | 8.8e-118 |
| std T3 row ('std',0) w4 | 5.0e-120 | 5.0e-120 | 5.0e-120 | 5.0e-120 |
| std T3 row ('std',1) w4 | 1.1e-120 | 2.0e-121 | 2.0e-121 | 1.1e-120 |
| shuffle h(a)h(b)=h(a,b)+h(b,a) (3,5) | 9.5e-212 | 1.8e-188 | 1.8e-188 | 9.6e-212 |
| shuffle (1,2) | 1.1e-188 | 1.7e-211 | 1.3e-211 | 1.1e-188 |
| random ℤ-combo of the 32 w3 words | 280.6 | 105.6 | 105.6 | 280.6 |

**Identified residual:** five-term at k=3 equals exactly **−4·π·i·log(1+√2)** (imag/(π·log(1+√2))
= −4 to 40 digits). √2 ↦ −√2 sends log(1+√2) to log(1−√2) = −log(1+√2) + π·i and the Rogers
correction terms pick up branch shifts — the relation is real-motivic, not form-level.
Conclusion: the trick holds for engine/form-level rows, FAILS for Rogers/R1–R5 at k=3,5.
→ Use single-embedding (k=1) LLL for the full 29-relation search; use 8-constraint LLL only to
identify the form-level subset (and report which are invariant).

## Step 2 — values for the 110 free words
- `vals110.pickle`: 110 words × k∈{1,3,5,7}, acb prec 1200 bits, cap 2400, N=2000, strings 350 digits.
- Convergence certificates per k: max|N2000−N1000| = 3.6e-190 (k=1), 1.5e-192 (k=3,5), 3.6e-190 (k=7);
  max|N2000−N1400| = 2.2e-264 (k=1), 8.3e-267 (k=3,5,7). Worst words: (7,7,7,0), (7,7,5,0).
  ⇒ values certified ≥ 258 digits, self-consistency ~369 digits.
- Independent check: nested mp.quad for G((7,4,7,7);1) matches engine to 1.62e-38 (596 s).

## Steps 3+ — pending
LLL search (single-embedding at ~230 digits; 8-constraint cross-check), D-target membership,
exact certificate, structure analysis. See step_NN.py and progress.log.
'''
open(WORK+'/result.md','w').write(md)
import datetime
with open(WORK+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} step16-33: setup evsetup; Galois trick tested = HALF-FAILS (Rogers/R1-5 fail k=3,5 by -4*pi*i*log(1+sqrt2) type terms); vals110.pickle (110 words x 4 embeddings, N=2000 cert 258d); quad check (7,4,7,7) 1.6e-38; LLL pipeline validated on w3 case; result.md early draft\n")
print('result.md written')
