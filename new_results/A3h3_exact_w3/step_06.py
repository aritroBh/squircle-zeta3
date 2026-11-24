import pickle
from pathlib import Path
from fractions import Fraction
import math
from flint import fmpz_mat, fmpz
from collections import defaultdict
import sys

def build_matrix_from_rows(word_to_idx, rows_data):
    """
    Build A and B matrices from rows_data where coefficients are Fractions or [x, y] pairs.
    Returns A, B where M = A + √2·B
    """
    n_words = len(word_to_idx)
    n_rows = len(rows_data)
    
    # Initialize matrices with zeros
    A = [[Fraction(0) for _ in range(n_words)] for _ in range(n_rows)]
    B = [[Fraction(0) for _ in range(n_words)] for _ in range(n_rows)]
    
    for i, row in enumerate(rows_data):
        if isinstance(row, list):
            # List of (word, coefficient) pairs
            for word, coeff in row:
                if word in word_to_idx:
                    j = word_to_idx[word]
                    if isinstance(coeff, Fraction):
                        A[i][j] += coeff
                    elif isinstance(coeff, (int, float)):
                        A[i][j] += Fraction(coeff)
                    elif isinstance(coeff, list) and len(coeff) == 2:
                        # [x, y] where x + √2·y
                        A[i][j] += Fraction(coeff[0])
                        B[i][j] += Fraction(coeff[1])
        elif isinstance(row, dict):
            # Dict word → [x, y] or (x, y)
            for word, coeff in row.items():
                if word in word_to_idx:
                    j = word_to_idx[word]
                    if isinstance(coeff, list) and len(coeff) == 2:
                        A[i][j] += Fraction(coeff[0])
                        B[i][j] += Fraction(coeff[1])
                    elif isinstance(coeff, tuple) and len(coeff) == 2:
                        A[i][j] += Fraction(coeff[0])
                        B[i][j] += Fraction(coeff[1])
        elif isinstance(row, tuple):
            # Tuple of (word, coefficient) pairs
            # Assuming it's a flat tuple: (word1, coeff1, word2, coeff2, ...)
            for k in range(0, len(row), 2):
                if k+1 < len(row):
                    word = row[k]
                    coeff = row[k+1]
                    if word in word_to_idx:
                        j = word_to_idx[word]
                        if isinstance(coeff, Fraction):
                            A[i][j] += coeff
                        elif isinstance(coeff, (int, float)):
                            A[i][j] += Fraction(coeff)
        else:
            print(f"Warning: row {i} has unexpected type {type(row)}")
    
    return A, B

def compute_qsqrt2_rank(A, B):
    """
    Compute rank over Q(√2) using the block matrix formula.
    A, B are lists of lists of Fractions.
    """
    n_rows = len(A)
    n_cols = len(A[0]) if n_rows > 0 else 0
    
    # Build block matrix [[A, 2B], [B, A]] with cleared denominators
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
            if d != 0:  # safety check
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
    
    # Rank over Q(√2) is rank/2 (should be integer)
    if rank % 2 != 0:
        print(f"Warning: block matrix rank {rank} is not even!")
    
    return rank // 2

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
print(f"Number of rows: {len(rows_data)}")
print(f"First row type: {type(rows_data[0])}")
print(f"First row length: {len(rows_data[0]) if isinstance(rows_data[0], (list, tuple)) else 'N/A'}")
print(f"First few elements of first row: {rows_data[0][:4] if isinstance(rows_data[0], (list, tuple)) else rows_data[0]}")

A1, B1 = build_matrix_from_rows(word_to_idx, rows_data)

print(f"\nA1 dimensions: {len(A1)} x {len(A1[0]) if A1 else 0}")
print(f"B1 dimensions: {len(B1)} x {len(B1[0]) if B1 else 0}")

# Check that B should be all zeros (since E1 has no √2 coefficients)
non_zero_B = 0
for i in range(len(B1)):
    for j in range(len(B1[0])):
        if B1[i][j] != 0:
            non_zero_B += 1
print(f"Non-zero B coefficients: {non_zero_B}")

# Compute rank
print("\nComputing rank of E1 over Q(√2)...")
rank_E1 = compute_qsqrt2_rank(A1, B1)
print(f"Exact rank of E1: {rank_E1}")
print(f"Expected rank: 542")
print(f"Match: {rank_E1 == 542}")