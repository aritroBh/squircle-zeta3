import pickle
from pathlib import Path
from fractions import Fraction
import math
from flint import fmpz_mat, fmpz
import datetime

# Reuse functions
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

def build_vector_from_relation(word_to_idx, relation_vector, freewords):
    """Build a full vector from relation_vector on freewords."""
    n_words = len(word_to_idx)
    vec = [Fraction(0) for _ in range(n_words)]
    
    for i, coeff in enumerate(relation_vector):
        word = freewords[i]
        if word in word_to_idx:
            j = word_to_idx[word]
            vec[j] += Fraction(coeff)
    
    return vec

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

base_path = Path('research/08_theory/')

print("Step 4: For each i, exact rank of E2 + {R_i}")
print("Expected rank: 545 (proves R_i is not in Q(√2)-span of E2)")

# Load data
with open(base_path / 'A3f3_level8_compute/weight3_with_lifts.pickle', 'rb') as f:
    w3_data = pickle.load(f)

with open(base_path / 'A3f3_level8_compute/k3_newHrows.pickle', 'rb') as f:
    k3_data = pickle.load(f)

with open(base_path / 'A3f4_glanois_basis/five_relations.pickle', 'rb') as f:
    rel_data = pickle.load(f)

# Build word to index mapping
words = w3_data['cols']
word_to_idx = {word: i for i, word in enumerate(words)}

# Get E1 rows and Kummer rows
E1_rows = w3_data['rows']
krows = k3_data['krows']

# Get relations
relations = rel_data['relations']
freewords = rel_data['freewords']

# E2 = E1 + Kummer rows
E2_rows = list(E1_rows) + list(krows)

# Convert each relation to the dict format needed by build_matrix_from_rows
def relation_to_row_dict(rel_vec):
    """Convert a relation vector to row dict format."""
    row_dict = {}
    vec = build_vector_from_relation(word_to_idx, rel_vec, freewords)
    for j, coeff in enumerate(vec):
        if coeff != 0:
            word = words[j]
            row_dict[word] = [coeff, Fraction(0)]  # [A, B] with B=0
    return row_dict

# Compute ranks for E2 + each R_i
results_step4 = []
for i, rel_vec in enumerate(relations):
    print(f"\nProcessing relation R_{i+1}...")
    
    # Create E2 + R_i
    current_rows = E2_rows + [relation_to_row_dict(rel_vec)]
    
    # Build matrices
    A, B = build_matrix_from_rows(word_to_idx, current_rows)
    
    # Compute rank
    rank = compute_qsqrt2_rank(A, B)
    results_step4.append(rank)
    
    print(f"  Rank of E2 + R_{i+1}: {rank}")
    print(f"  Expected: 545")
    print(f"  Rank increases by {rank - 544} (should be +1)")
    
    if rank == 545:
        print(f"  ✓ R_{i+1} is NOT in the Q(√2)-span of E2 (PROVED)")
    else:
        print(f"  ✗ Unexpected result")

# Log results
with open('progress.log', 'a') as log:
    log.write(f"{datetime.datetime.now()}: Step 4 completed\n")
    for i, rank in enumerate(results_step4):
        log.write(f"  R_{i+1}: rank = {rank} (expected 545), {'✓ PROVED' if rank == 545 else '✗ FAILED'}\n")

print("\n" + "="*60)
print("Step 4 Summary:")
for i, rank in enumerate(results_step4):
    status = "✓ PROVED" if rank == 545 else "✗ FAILED"
    print(f"R_{i+1}: rank = {rank}, expected = 545, {status}")