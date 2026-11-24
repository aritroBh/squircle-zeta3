# Create step01_rank_E1.py
step01_content = '''#!/usr/bin/env python3
"""Step 1: Exact rank of E1 = weight3_with_lifts rows (expect 542)"""

import pickle
from pathlib import Path
from fractions import Fraction
import math
from flint import fmpz_mat, fmpz

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
    
    return rank // 2

def main():
    base_path = Path('research/08_theory/')
    
    print("Loading weight3_with_lifts.pickle...")
    with open(base_path / 'A3f3_level8_compute/weight3_with_lifts.pickle', 'rb') as f:
        w3_data = pickle.load(f)
    
    # Build word to index mapping
    words = w3_data['cols']
    word_to_idx = {word: i for i, word in enumerate(words)}
    
    # Build matrices for E1
    print(f"Building matrices from {len(w3_data['rows'])} rows...")
    A1, B1 = build_matrix_from_rows(word_to_idx, w3_data['rows'])
    
    print(f"Matrix dimensions: {len(A1)} × {len(A1[0])}")
    
    # Compute rank
    print("Computing rank over Q(√2)...")
    rank_E1 = compute_qsqrt2_rank(A1, B1)
    
    print(f"\\nResult:")
    print(f"  Exact rank of E1: {rank_E1}")
    print(f"  Expected rank: 542")
    print(f"  Match: {rank_E1 == 542}")
    
    if rank_E1 == 542:
        print("\\n✓ Step 1: PROVED")
    else:
        print("\\n✗ Step 1: FAILED")

if __name__ == "__main__":
    main()
'''

with open('step01_rank_E1.py', 'w') as f:
    f.write(step01_content)
print("Created step01_rank_E1.py")