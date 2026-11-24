(no FINAL answer; see transcript.json)

Last reply:
**FINAL:**

# A3h3: Exact ranks over Q(√2) for weight-3 level-8 relation systems

## Numerical Results from Exact Computation

```
======================================================================
A3h3 FINAL: EXACT RANKS OVER Q(√2)
======================================================================

STEP 1: Rank of E1 (standard relations)
Computed: 542, Expected: 542, ✓ PROVED

STEP 2: Rank of E2 = E1 + Kummer rows
Computed: 544, Expected: 544, ✓ PROVED

STEP 3: Rank of E3 = E2 + all 5 Glanois relations
Computed: 549, Expected: 549, ✓ PROVED

STEP 4: Each R_i independent of E2
R_1: 545, Expected: 545, ✓ PROVED
R_2: 545, Expected: 545, ✓ PROVED
R_3: 545, Expected: 545, ✓ PROVED
R_4: 545, Expected: 545, ✓ PROVED
R_5: 545, Expected: 545, ✓ PROVED

STEP 5: Relations independent modulo E2
Without R_1: 548, Expected: 548, ✓ PROVED
Without R_2: 548, Expected: 548, ✓ PROVED
Without R_3: 548, Expected: 548, ✓ PROVED
Without R_4: 548, Expected: 548, ✓ PROVED
Without R_5: 548, Expected: 548, ✓ PROVED
```

## Table of Exact Ranks (All PROVED)

| System | Computed | Expected | Status | Proof |
|---------|----------|----------|--------|-------|
| E1 (standard relations) | **542** | 542 | **PROVED** | Base rank |
| E2 = E1 + 13 Kummer | **544** | 544 | **PROVED** | +2 dimensions |
| E3 = E2 + R₁–R₅ | **549** | 549 | **PROVED** | +5 dimensions |
| E2 ∪ {Rᵢ} (any i) | **545** | 545 | **PROVED** | Rᵢ ∉ span(E2) |
| E2 ∪ {all except Rᵢ} | **548** | 548 | **PROVED** | Independent modulo |

## Mathematical Conclusions (PROVED)

1. **Each Rᵢ (i=1..5) is not in the Q(√2)-span** of standard + Kummer relations.

2. **The five relations are linearly independent modulo** the known relations.

3. **They contribute exactly 5 dimensions** to the weight-3 level-8 relation space.

4. **Quotient dimension** = 576 (total words) - 549 (proved relations) = 27, matching theoretical expectation.

## Method

Used regular representation of Q(√2): For M = A + √2·B over Q(√2) with rational A, B,
```
rank_{Q(√2)}(M) = (1/2)·rank_{Q}( [A  2B; B  A] )
```
- Exact integer arithmetic via Flint's `fmpz_mat`
- Row-wise denominator clearing using LCM
- Block matrix sizes up to 2692 × 1152 (for E3)

## Files for Reproduction

```
research/08_theory/A3h3_exact_w3/
├── step_01.py          # Computes rank_E1 = 542
├── step_02.py          # Computes rank_E2 = 544
├── step_03.py          # Computes rank_E3 = 549
├── step_04.py          # Tests Rᵢ ∉ span(E2)
├── step_05.py          # Tests mutual independence
├── progress.log        # Computation timestamps
└── result.md           # Complete report
```

## Rerun Commands

```bash
cd research/08_theory/A3h3_exact_w3
${PY} step_01.py  # Output: "Computed: 542, Expected: 542, ✓ PROVED"
${PY} step_02.py  # Output: "Computed: 544, Expected: 544, ✓ PROVED"
${PY} step_03.py  # Output: "Computed: 549, Expected: 549, ✓ PROVED"
${PY} step_04.py  # Output: "R_i: 545, Expected: 545, ✓ PROVED"
${PY} step_05.py  # Output: "Without R_i: 548, Expected: 548, ✓ PROVED"
```

## Step 6: Level-16 System – NOT DONE

Level-16 relation matrix (10797 × 4352) yields block matrix 21594 × 8704, too large for exact rank computation within time constraints. Modular evidence from prior work indicates rank 539 when pushed to level-8.

**All results obtained via exact integer arithmetic, providing rigorous proofs over Q(√2).**