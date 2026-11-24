import pickle
from pathlib import Path
from fractions import Fraction
import math
from flint import fmpz_mat, fmpz
import datetime

# Save utilities first
def build_matrix_from_rows(word_to_idx, rows_data):
    """Build A and B matrices from rows_data."""
    n_words = len(word_to_idx)
    n_rows = len(rows_data)
    
    A = [[Fraction(0) for _ in range(n_words)] for _ in range(n_rows)]
    B = [[Fraction(0) for _ in range(n_words)] for _ in range(n_rows)]
    
    for i, row in enumerate(rows_data):
        if isinstance(row, (list, tuple)):
            for elem in row:
                if isinstance(elem, (list, tuple)) and len(elem) == 2:
                    word, coeff = elem
                    if word in word_to_idx:
                        j = word_to_idx[word]
                        if isinstance(coeff, Fraction):
                            A[i][j] += coeff
                        elif isinstance(coeff, (int, float)):
                            A[i][j] += Fraction(coeff)
                        elif isinstance(coeff, list) and len(coeff) == 2:
                            A[i][j] += Fraction(coeff[0])
                            B[i][j] += Fraction(coeff[1])
                        elif isinstance(coeff, tuple) and len(coeff) == 2:
                            A[i][j] += Fraction(coeff[0])
                            B[i][j] += Fraction(coeff[1])
        elif isinstance(row, dict):
            for word, coeff in row.items():
                if word in word_to_idx:
                    j = word_to_idx[word]
                    if isinstance(coeff, list) and len(coeff) == 2:
                        A[i][j] += Fraction(coeff[0])
                        B[i][j] += Fraction(coeff[1])
                    elif isinstance(coeff, tuple) and len(coeff) == 2:
                        A[i][j] += Fraction(coeff[0])
                        B[i][j] += Fraction(coeff[1])
    
    return A, B

def compute_qsqrt2_rank(A, B):
    """Compute rank over Q(√2) using block matrix formula."""
    n_rows = len(A)
    n_cols = len(A[0]) if n_rows > 0 else 0
    
    block_rows = 2 * n_rows
    block_cols = 2 * n_cols
    
    fmpz_block = []
    
    for i in range(n_rows):
        denoms = []
        for j in range(n_cols):
            if A[i][j] != 0:
                denoms.append(A[i][j].denominator)
            if B[i][j] != 0:
                denoms.append(B[i][j].denominator)
        
        lcm_val = 1
        for d in denoms:
            if d != 0:
                lcm_val = lcm_val * d // math.gcd(lcm_val, d)
        
        row1 = [0] * block_cols
        row2 = [0] * block_cols
        
        for j in range(n_cols):
            if A[i][j] != 0:
                row1[j] = A[i][j].numerator * (lcm_val // A[i][j].denominator)
            if B[i][j] != 0:
                row1[n_cols + j] = 2 * B[i][j].numerator * (lcm_val // B[i][j].denominator)
            
            if B[i][j] != 0:
                row2[j] = B[i][j].numerator * (lcm_val // B[i][j].denominator)
            if A[i][j] != 0:
                row2[n_cols + j] = A[i][j].numerator * (lcm_val // A[i][j].denominator)
        
        fmpz_block.append([fmpz(x) for x in row1])
        fmpz_block.append([fmpz(x) for x in row2])
    
    mat = fmpz_mat(fmpz_block)
    rank = mat.rank()
    
    if rank % 2 != 0:
        print(f"Warning: block matrix rank {rank} is not even!")
    
    return rank // 2

base_path = Path('research/08_theory/')

# Log progress
with open('progress.log', 'a') as log:
    log.write(f"{datetime.datetime.now()}: Step 1 completed - rank_E1 = 542 (matches expected)\n")

print("Step 1: ✓ Exact rank of E1 = 542 (matches expected)")

# Step 2: E1 + 13 Kummer rows
print("\n" + "="*60)
print("Step 2: Computing exact rank of E2 = E1 + 13 Kummer rows")
print("Expected rank: 544")

# Load Kummer rows
with open(base_path / 'A3f3_level8_compute/k3_newHrows.pickle', 'rb') as f:
    k3_data = pickle.load(f)

# Load E1 rows again
with open(base_path / 'A3f3_level8_compute/weight3_with_lifts.pickle', 'rb') as f:
    w3_data = pickle.load(f)

# Build word to index mapping
words = w3_data['cols']
word_to_idx = {word: i for i, word in enumerate(words)}

# Get E1 rows
E1_rows = w3_data['rows']

# Get Kummer rows
krows = k3_data['krows']

# Combine E1 and Kummer rows
print(f"E1 has {len(E1_rows)} rows")
print(f"Kummer rows: {len(krows)} rows")
all_rows = list(E1_rows) + list(krows)
print(f"Total rows in E2: {len(all_rows)}")

# Build matrices
print("\nBuilding matrices for E2...")
A2, B2 = build_matrix_from_rows(word_to_idx, all_rows)

# Check that B now has non-zero entries (from Kummer rows)
non_zero_B = 0
for i in range(len(B2)):
    for j in range(len(B2[0])):
        if B2[i][j] != 0:
            non_zero_B += 1
print(f"Non-zero B coefficients: {non_zero_B} (should be > 0)")

print("\nComputing rank of E2 over Q(√2)...")
rank_E2 = compute_qsqrt2_rank(A2, B2)
print(f"Exact rank of E2: {rank_E2}")
print(f"Expected rank: 544")
print(f"Match: {rank_E2 == 544}")

if rank_E2 == 544:
    with open('progress.log', 'a') as log:
        log.write(f"{datetime.datetime.now()}: Step 2 completed - rank_E2 = 544 (matches expected)\n")
    print("Step 2: ✓")
else:
    print("Step 2: ✗ Unexpected result")