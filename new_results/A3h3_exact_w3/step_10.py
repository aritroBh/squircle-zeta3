import pickle
from pathlib import Path
from fractions import Fraction
import math
from flint import fmpz_mat, fmpz
import datetime

# Reuse functions from previous steps
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
    
    if rank % 2 != 0:
        print(f"Warning: block matrix rank {rank} is not even!")
    
    return rank // 2

base_path = Path('research/08_theory/')

print("Step 3: Computing exact rank of E3 = E2 + R1..R5")
print("Expected rank: 549 (576 - 27)")

# Load all data
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

# Build relation vectors
relations = rel_data['relations']
freewords = rel_data['freewords']

print(f"\nRelations: {len(relations)}")
print(f"Freewords: {len(freewords)}")
print(f"First freeword: {freewords[0]}")
print(f"First relation first 5 coefficients: {relations[0][:5]}")

# Convert relations to A, B format (all A, no √2 component)
rel_vectors = []
for rel_vec in relations:
    vec = build_vector_from_relation(word_to_idx, rel_vec, freewords)
    rel_vectors.append(vec)

print(f"\nBuilt {len(rel_vectors)} relation vectors")

# Combine all rows: E1 + Kummer + relations
all_rows = list(E1_rows) + list(krows) + rel_vectors

# But we need to convert relation vectors to the same format as other rows
# The build_matrix_from_rows expects (word, coeff) pairs or dicts
# So we need to create a list of dicts for the relation vectors
rel_rows = []
for vec in rel_vectors:
    row_dict = {}
    for j, coeff in enumerate(vec):
        if coeff != 0:
            word = words[j]
            row_dict[word] = [coeff, Fraction(0)]  # [A, B] with B=0
    rel_rows.append(row_dict)

# Now combine everything
all_rows_combined = list(E1_rows) + list(krows) + rel_rows
print(f"\nTotal rows in E3: {len(all_rows_combined)}")

# Build matrices
print("\nBuilding matrices for E3...")
A3, B3 = build_matrix_from_rows(word_to_idx, all_rows_combined)

print(f"A3 dimensions: {len(A3)} x {len(A3[0])}")
print(f"B3 dimensions: {len(B3)} x {len(B3[0])}")

# Compute rank
print("\nComputing rank of E3 over Q(√2)...")
rank_E3 = compute_qsqrt2_rank(A3, B3)
print(f"Exact rank of E3: {rank_E3}")
print(f"Expected rank: 549")
print(f"Difference from E2: {rank_E3 - 544} (should be +5)")

if rank_E3 == 549:
    with open('progress.log', 'a') as log:
        log.write(f"{datetime.datetime.now()}: Step 3 completed - rank_E3 = 549 (matches expected)\n")
    print("Step 3: ✓")
else:
    print("Step 3: ✗ Unexpected result")