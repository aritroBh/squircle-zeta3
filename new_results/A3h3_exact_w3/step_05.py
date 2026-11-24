import pickle
from pathlib import Path
from fractions import Fraction
import math
from flint import fmpz_mat, fmpz
from collections import defaultdict

base_path = Path('research/08_theory/')

# Load data
print("Loading data...")
with open(base_path / 'A3f3_level8_compute/weight3_with_lifts.pickle', 'rb') as f:
    w3_data = pickle.load(f)

# Build word to index mapping
words = w3_data['cols']  # 576 words
word_to_idx = {word: i for i, word in enumerate(words)}
print(f"Word dictionary: {len(word_to_idx)} words")

# Build A, B matrices for E1
print("\nBuilding matrices for E1 (weight3_with_lifts)...")
rows_data = w3_data['rows']
A1, B1 = build_matrix_from_rows(word_to_idx, rows_data)

print(f"A1 dimensions: {len(A1)} x {len(A1[0]) if A1 else 0}")
print(f"B1 dimensions: {len(B1)} x {len(B1[0]) if B1 else 0}")

# Compute rank
print("\nComputing rank of E1 over Q(√2)...")
rank_E1 = compute_qsqrt2_rank(A1, B1)
print(f"Exact rank of E1: {rank_E1}")
print(f"Expected rank: 542")
print(f"Match: {rank_E1 == 542}")

# Save this function for later use
def compute_qsqrt2_rank_optimized(A, B):
    """Optimized version of the rank computation."""
    n_rows = len(A)
    n_cols = len(A[0]) if n_rows > 0 else 0
    
    # Build block matrix directly with cleared denominators
    block_rows = 2 * n_rows
    block_cols = 2 * n_cols
    
    # We'll build the matrix in flint format directly
    fmpz_block = []
    
    for i in range(n_rows):
        # Find LCM for this row
        denoms = []
        for j in range(n_cols):
            if A[i][j] != 0:
                denoms.append(A[i][j].denominator)
            if B[i][j] != 0:
                denoms.append(B[i][j].denominator)
        
        lcm_val = 1
        for d in denoms:
            lcm_val = lcm_val * d // math.gcd(lcm_val, d)
        
        # Build two rows of the block matrix
        row1 = [0] * block_cols
        row2 = [0] * block_cols
        
        for j in range(n_cols):
            # Top row: A and 2B
            if A[i][j] != 0:
                row1[j] = A[i][j].numerator * (lcm_val // A[i][j].denominator)
            if B[i][j] != 0:
                row1[n_cols + j] = 2 * B[i][j].numerator * (lcm_val // B[i][j].denominator)
            
            # Bottom row: B and A
            if B[i][j] != 0:
                row2[j] = B[i][j].numerator * (lcm_val // B[i][j].denominator)
            if A[i][j] != 0:
                row2[n_cols + j] = A[i][j].numerator * (lcm_val // A[i][j].denominator)
        
        fmpz_block.append([fmpz(x) for x in row1])
        fmpz_block.append([fmpz(x) for x in row2])
    
    mat = fmpz_mat(fmpz_block)
    rank = mat.rank()
    return rank // 2

# Test with optimized function
print("\nComputing rank with optimized function...")
rank_E1_opt = compute_qsqrt2_rank_optimized(A1, B1)
print(f"Exact rank of E1 (optimized): {rank_E1_opt}")
print(f"Match with previous: {rank_E1 == rank_E1_opt}")