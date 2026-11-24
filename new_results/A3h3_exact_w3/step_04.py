from fractions import Fraction
from collections import defaultdict
import math
from flint import fmpz_mat, fmpz

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
    
    return A, B

def build_vector_from_relation(word_to_idx, relation_vector, freewords):
    """
    Build a vector for a relation given relation_vector (list of coefficients)
    and freewords (list of words).
    """
    n_words = len(word_to_idx)
    vec = [Fraction(0) for _ in range(n_words)]
    
    for i, coeff in enumerate(relation_vector):
        word = freewords[i]
        if word in word_to_idx:
            j = word_to_idx[word]
            vec[j] += Fraction(coeff)
    
    return vec

def compute_qsqrt2_rank(A, B):
    """
    Compute rank over Q(√2) using the block matrix formula.
    A, B are lists of lists of Fractions.
    """
    n_rows = len(A)
    n_cols = len(A[0]) if n_rows > 0 else 0
    
    # Clear denominators row by row
    # Find LCM for each row
    row_lcms = []
    for i in range(n_rows):
        denoms = []
        for j in range(n_cols):
            if A[i][j] != 0:
                denoms.append(A[i][j].denominator)
            if B[i][j] != 0:
                denoms.append(B[i][j].denominator)
        
        if denoms:
            # Compute LCM of denominators
            lcm_val = 1
            for d in denoms:
                lcm_val = lcm_val * d // math.gcd(lcm_val, d)
            row_lcms.append(lcm_val)
        else:
            row_lcms.append(1)
    
    # Build block matrix [[A, 2B], [B, A]] with cleared denominators
    block_rows = 2 * n_rows
    block_cols = 2 * n_cols
    
    # Initialize integer block matrix
    block_matrix = [[0 for _ in range(block_cols)] for _ in range(block_rows)]
    
    for i in range(n_rows):
        lcm = row_lcms[i]
        # Top block: A(i) and 2B(i)
        for j in range(n_cols):
            # A part
            block_matrix[i][j] = int(A[i][j].numerator * (lcm // A[i][j].denominator) if A[i][j] != 0 else 0)
            # 2B part
            b_val = B[i][j]
            block_matrix[i][n_cols + j] = int(2 * b_val.numerator * (lcm // b_val.denominator) if b_val != 0 else 0)
        
        # Bottom block: B(i) and A(i)
        for j in range(n_cols):
            # B part
            b_val = B[i][j]
            block_matrix[n_rows + i][j] = int(b_val.numerator * (lcm // b_val.denominator) if b_val != 0 else 0)
            # A part
            block_matrix[n_rows + i][n_cols + j] = int(A[i][j].numerator * (lcm // A[i][j].denominator) if A[i][j] != 0 else 0)
    
    # Convert to flint matrix and compute rank
    fmpz_block = []
    for row in block_matrix:
        fmpz_row = [fmpz(x) for x in row]
        fmpz_block.append(fmpz_row)
    
    mat = fmpz_mat(fmpz_block)
    rank = mat.rank()
    
    # Rank over Q(√2) is rank/2 (should be integer)
    if rank % 2 != 0:
        print(f"Warning: block matrix rank {rank} is not even!")
    
    return rank // 2

# Test with a small example
print("Testing utility functions...")
word_to_idx_test = {(1, 2): 0, (2, 1): 1, (3, 4): 2}
rows_test = [
    [((1, 2), Fraction(1, 2)), ((2, 1), Fraction(1, 3))],
    [((3, 4), Fraction(2, 3))]
]

A_test, B_test = build_matrix_from_rows(word_to_idx_test, rows_test)
print(f"A_test[0]: {A_test[0]}")
print(f"A_test[1]: {A_test[1]}")
print(f"B_test[0]: {B_test[0]}")
print(f"B_test[1]: {B_test[1]}")

# Test with Kummer-style row
krow_test = [{(1, 2): [Fraction(1, 2), Fraction(1, 3)], (2, 1): [Fraction(1, 4), Fraction(1, 5)]}]
A_k, B_k = build_matrix_from_rows(word_to_idx_test, krow_test)
print(f"\nKummer test:")
print(f"A_k[0]: {A_k[0]}")
print(f"B_k[0]: {B_k[0]}")
