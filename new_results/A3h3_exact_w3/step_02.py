import sympy as sp
import math
from flint import fmpz_mat, fmpq_mat, fmpz, fmpq

# Test the rank formula for Q(√2)
print("Testing rank formula for Q(√2)")
print("=" * 50)

# Simple test matrix over Q(√2)
# M = [[1, √2], [√2, 2]] which has rank 1 since second row = √2 × first row
# Represent as A + √2·B where A,B are rational
M = [[1, sp.sqrt(2)], [sp.sqrt(2), 2]]

# Convert to A, B matrices
A = [[1, 0], [0, 2]]  # rational part
B = [[0, 1], [1, 0]]  # √2 coefficient part

print("Matrix M (over Q(√2)):")
print(M)
print("\nA (rational part):")
for row in A:
    print(row)
print("\nB (√2 coefficient part):")
for row in B:
    print(row)

# Build the block matrix [[A, 2B], [B, A]]
rows = len(A)
cols = len(A[0])
block_rows = 2*rows
block_cols = 2*cols

block_matrix = []
for i in range(rows):
    row = A[i] + [2*x for x in B[i]]
    block_matrix.append(row)
for i in range(rows):
    row = B[i] + A[i]
    block_matrix.append(row)

print("\nBlock matrix [[A, 2B], [B, A]]:")
for row in block_matrix:
    print(row)

# Convert to flint matrix and compute rank
fmpz_block = []
for row in block_matrix:
    fmpz_row = []
    for x in row:
        fmpz_row.append(fmpz(int(x)))
    fmpz_block.append(fmpz_row)

mat = fmpz_mat(fmpz_block)
rank = mat.rank()
print(f"\nRank of block matrix: {rank}")
print(f"Rank over Q(√2) = rank/2 = {rank/2}")

# Let's also test with a full rank 2x2 matrix
print("\n" + "=" * 50)
print("Testing with full rank 2x2 matrix:")
# M = [[1, √2], [1, 1]]
A2 = [[1, 0], [1, 1]]
B2 = [[0, 1], [0, 0]]

block_matrix2 = []
for i in range(2):
    row = A2[i] + [2*x for x in B2[i]]
    block_matrix2.append(row)
for i in range(2):
    row = B2[i] + A2[i]
    block_matrix2.append(row)

print("Block matrix:")
for row in block_matrix2:
    print(row)

fmpz_block2 = []
for row in block_matrix2:
    fmpz_row = []
    for x in row:
        fmpz_row.append(fmpz(int(x)))
    fmpz_block2.append(fmpz_row)

mat2 = fmpz_mat(fmpz_block2)
rank2 = mat2.rank()
print(f"Rank of block matrix: {rank2}")
print(f"Rank over Q(√2) = rank/2 = {rank2/2} (should be 2)")