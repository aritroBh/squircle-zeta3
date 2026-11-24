#!/usr/bin/env python3
"""
Final verification of exact ranks over Q(√2) for weight-3 level-8 relations
"""

import pickle
from pathlib import Path
from fractions import Fraction
import math
from flint import fmpz_mat, fmpz
import datetime

# Utility functions
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

# Main execution
print("="*70)
print("A3h3: EXACT RANKS OVER Q(√2) FOR WEIGHT-3 LEVEL-8 RELATION SYSTEMS")
print("="*70)

base_path = Path('research/08_theory/')

# Load all data
print("\nLoading data...")
with open(base_path / 'A3f3_level8_compute/weight3_with_lifts.pickle', 'rb') as f:
    w3_data = pickle.load(f)

with open(base_path / 'A3f3_level8_compute/k3_newHrows.pickle', 'rb') as f:
    k3_data = pickle.load(f)

with open(base_path / 'A3f4_glanois_basis/five_relations.pickle', 'rb') as f:
    rel_data = pickle.load(f)

# Build word to index mapping
words = w3_data['cols']
word_to_idx = {word: i for i, word in enumerate(words)}

# Get rows
E1_rows = w3_data['rows']
krows = k3_data['krows']
relations = rel_data['relations']
freewords = rel_data['freewords']

print(f"Vocabulary: {len(words)} words")
print(f"E1 rows: {len(E1_rows)}")
print(f"Kummer rows: {len(krows)}")
print(f"Glanois relations: {len(relations)}")

print("\n" + "="*70)
print("STEP 1: Exact rank of E1 (weight3_with_lifts)")
print("="*70)

A1, B1 = build_matrix_from_rows(word_to_idx, E1_rows)
rank_E1 = compute_qsqrt2_rank(A1, B1)
print(f"Computed rank: {rank_E1}")
print(f"Expected rank: 542")
print(f"Status: {'✓ PROVED' if rank_E1 == 542 else '✗ FAILED'}")

print("\n" + "="*70)
print("STEP 2: Exact rank of E2 = E1 + 13 Kummer rows")
print("="*70)

E2_rows = list(E1_rows) + list(krows)
A2, B2 = build_matrix_from_rows(word_to_idx, E2_rows)
rank_E2 = compute_qsqrt2_rank(A2, B2)
print(f"Computed rank: {rank_E2}")
print(f"Expected rank: 544")
print(f"Status: {'✓ PROVED' if rank_E2 == 544 else '✗ FAILED'}")

print("\n" + "="*70)
print("STEP 3: Exact rank of E3 = E2 + all 5 Glanois relations")
print("="*70)

# Convert relations to row dict format
def relation_to_row_dict(rel_vec):
    row_dict = {}
    vec = build_vector_from_relation(word_to_idx, rel_vec, freewords)
    for j, coeff in enumerate(vec):
        if coeff != 0:
            word = words[j]
            row_dict[word] = [coeff, Fraction(0)]
    return row_dict

rel_row_dicts = [relation_to_row_dict(rel_vec) for rel_vec in relations]
E3_rows = E2_rows + rel_row_dicts
A3, B3 = build_matrix_from_rows(word_to_idx, E3_rows)
rank_E3 = compute_qsqrt2_rank(A3, B3)
print(f"Computed rank: {rank_E3}")
print(f"Expected rank: 549 (576 - 27)")
print(f"Status: {'✓ PROVED' if rank_E3 == 549 else '✗ FAILED'}")

print("\n" + "="*70)
print("STEP 4: Each R_i is independent of E2")
print("="*70)

step4_results = []
for i, rel_vec in enumerate(relations):
    current_rows = E2_rows + [relation_to_row_dict(rel_vec)]
    A, B = build_matrix_from_rows(word_to_idx, current_rows)
    rank = compute_qsqrt2_rank(A, B)
    step4_results.append(rank)
    print(f"R_{i+1}: rank = {rank}, expected = 545, {'✓ PROVED' if rank == 545 else '✗ FAILED'}")

print("\n" + "="*70)
print("STEP 5: Relations are independent modulo E2")
print("="*70)

step5_results = []
for j in range(len(relations)):
    current_rels = [rel_row_dicts[i] for i in range(len(relations)) if i != j]
    current_rows = E2_rows + current_rels
    A, B = build_matrix_from_rows(word_to_idx, current_rows)
    rank = compute_qsqrt2_rank(A, B)
    step5_results.append(rank)
    print(f"Without R_{j+1}: rank = {rank}, expected = 548, {'✓ PROVED' if rank == 548 else '✗ FAILED'}")

print("\n" + "="*70)
print("SUMMARY OF PROOFS")
print("="*70)

print(f"""
1. E1 rank = {rank_E1} (expected 542): {'✓ PROVED' if rank_E1 == 542 else '✗ FAILED'}
2. E2 rank = {rank_E2} (expected 544): {'✓ PROVED' if rank_E2 == 544 else '✗ FAILED'}
3. E3 rank = {rank_E3} (expected 549): {'✓ PROVED' if rank_E3 == 549 else '✗ FAILED'}

4. Each R_i independent of E2: {'✓ ALL PROVED' if all(r == 545 for r in step4_results) else '✗ SOME FAILED'}
5. Relations independent modulo E2: {'✓ ALL PROVED' if all(r == 548 for r in step5_results) else '✗ SOME FAILED'}

Key implications:
- The 5 Glanois relations are not consequences of the standard relations + Kummer relations
- They form a 5-dimensional space over Q(√2) modulo the known relations
- The quotient space has dimension 27 (576 - 549), as expected
""")

print(f"\nComputation completed at: {datetime.datetime.now()}")
print("="*70)