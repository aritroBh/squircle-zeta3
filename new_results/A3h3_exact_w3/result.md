# A3h3: Exact ranks over Q(√2) for weight-3 level-8 relation systems

## Summary of Results (PROVED)

All computations were performed exactly using the regular representation 
of Q(√2): for a matrix M = A + √2·B over Q(√2), the rank over Q(√2) equals 
(rank over Q of the block matrix [[A, 2B], [B, A]]) divided by 2.

### Step 1: Standard relations E1 (`weight3_with_lifts`)
- **Computed rank**: 542
- **Expected rank**: 542
- **Status**: ✓ PROVED

### Step 2: E1 + 13 Kummer relations (E2)
- **Computed rank**: 544
- **Expected rank**: 544
- **Status**: ✓ PROVED
- **Interpretation**: The 13 Kummer relations add exactly 2 dimensions to the space

### Step 3: E2 + all 5 Glanois relations (E3)
- **Computed rank**: 549
- **Expected rank**: 549 (576 - 27)
- **Status**: ✓ PROVED
- **Interpretation**: The 5 Glanois relations add exactly 5 dimensions to the space

### Step 4: Each R_i is independent of E2
For each i = 1..5, rank(E2 ∪ {R_i}) was computed:
- R₁: rank = 545 (expected 545): ✓ PROVED (R_1 is NOT in the Q(√2)-span of E2)
- R₂: rank = 545 (expected 545): ✓ PROVED (R_2 is NOT in the Q(√2)-span of E2)
- R₃: rank = 545 (expected 545): ✓ PROVED (R_3 is NOT in the Q(√2)-span of E2)
- R₄: rank = 545 (expected 545): ✓ PROVED (R_4 is NOT in the Q(√2)-span of E2)
- R₅: rank = 545 (expected 545): ✓ PROVED (R_5 is NOT in the Q(√2)-span of E2)

**Conclusion**: Each of the five Glanois relations is linearly independent 
from the space spanned by the standard and Kummer relations over Q(√2).

### Step 5: The five relations are independent modulo E2
For each i = 1..5, rank(E2 ∪ {all R_j except R_i}) was computed:
- Without R₁: rank = 548 (expected 548): ✓ PROVED
- Without R₂: rank = 548 (expected 548): ✓ PROVED
- Without R₃: rank = 548 (expected 548): ✓ PROVED
- Without R₄: rank = 548 (expected 548): ✓ PROVED
- Without R₅: rank = 548 (expected 548): ✓ PROVED

**Conclusion**: All five Glanois relations are linearly independent modulo E2.
Together with Step 3 (rank = 549), this proves that the five relations form 
a 5-dimensional space modulo the standard and Kummer relations.

## Method Details

### Regular representation of Q(√2)
For an element a + b√2 ∈ Q(√2), its action on Q² is given by the matrix:
```
[ a  2b ]
[ b   a ]
```

For a matrix M = A + √2·B over Q(√2) with rational matrices A, B, 
we construct the block matrix:
```
[ A   2B ]
[ B    A ]
```

The rank of M over Q(√2) equals half the rank of this block matrix over Q.

### Implementation
- Coefficients were handled exactly using Python's `Fraction` type
- Denominators were cleared row-by-row using LCM
- Exact rank computation used Flint's `fmpz_mat` with integer arithmetic
- All computations were performed in exact integer arithmetic

### Matrix Sizes
- E1: 1328 × 576 (block: 2656 × 1152)
- E2: 1341 × 576 (block: 2682 × 1152)  
- E3: 1346 × 576 (block: 2692 × 1152)

## Files Produced

1. `progress.log` - Timestamped completion log of each step
2. All Python scripts used for computations

## Rerun Commands

The computations can be reproduced with:

```bash
cd research/08_theory/A3h3_exact_w3
```

### Step 1: Standard relations rank
```bash
${PY} step01_rank_E1.py
```

### Step 2: With Kummer relations  
```bash
${PY} step02_rank_E2.py
```

### Step 3: With all Glanois relations
```bash
${PY} step03_rank_E3.py
```

### Step 4: Individual independence
```bash
${PY} step04_individual_independence.py
```

### Step 5: Mutual independence
```bash
${PY} step05_mutual_independence.py
```

(Note: `${PY}` refers to `<project>/.venv/bin/python`)

## Step 6 (Optional): Level-16 system - NOT DONE

The level-16 system (`A3f6_level16/l16_w3.pickle`) has:
- Columns: 4352 level-16 words
- Rows: 10797 level-16 relations
- Block matrix size would be: 21594 × 8704

While the embedding of level-8 relations into level-16 has been verified 
(letter j ↦ 2j, with -1 fixed), computing the exact rank of such a large 
system was deemed infeasible within the time constraints (25 minutes).

Based on the results in `A3f6_level16/result.md`, it is known from modular 
computations that pushing level-16 relations down to level-8 gives rank 539, 
and that R₁–R₅ are not provable from level-16 standard relations.

**Status**: NOT DONE (matrix too large for exact rank computation within time limit)

---

*Report generated: 2026-10-22 14:10:26*
*All rank computations were performed using exact integer arithmetic via Flint's fmpz_mat.*

## Planner verification (W4 side session, 2026-10-22)
- DeepSeek V3.2 hit its turn limit (NO-FINAL) after writing this report. All five required steps finished; step 6 was
  not done.
- The planner recomputed steps 2–4 with independent code (`planner/planner_check_exact.py`). It uses the same regular
  representation [[A, 2B], [B, A]] and flint fmpz_mat.rank, and tests the block formula on a toy example first. The
  R1–R5 word order is read from five_relations.pickle['freewords'].
  - exact rank E2 = 544;
  - exact rank E2 ∪ {R_i} = 545 for each i;
  - exact rank E2 ∪ {R1..R5} = 549.
  This agrees with the table above.
- Conclusion (PROVED): R1–R5 are not in the ℚ(√2)-span of the weight-3 standard + lift + Kummer rows, and they are
  independent modulo that span. This makes the weight-3 negative statements exact rather than two-prime evidence.
