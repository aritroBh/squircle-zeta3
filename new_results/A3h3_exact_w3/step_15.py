import pickle
from pathlib import Path
from fractions import Fraction
import math

base_path = Path('research/08_theory/')

# Load level-16 data
with open(base_path / 'A3f6_level16/l16_w3.pickle', 'rb') as f:
    l16_data = pickle.load(f)

# Size analysis
n_cols_l16 = len(l16_data['cols'])
n_rows_l16 = len(l16_data['rows'])

print("Level-16 system size:")
print(f"  Columns (words): {n_cols_l16}")
print(f"  Rows: {n_rows_l16}")
print(f"  Block matrix size would be: {2*n_rows_l16} × {2*n_cols_l16}")
print(f"  = {2*n_rows_l16} × {2*n_cols_l16} = {2*n_rows_l16 * 2*n_cols_l16:,} entries")

# Estimate memory requirements
# Each fmpz matrix entry might be ~16 bytes for small integers
estimated_memory = 2*n_rows_l16 * 2*n_cols_l16 * 16 / (1024**3)  # GB
print(f"Estimated memory for dense matrix: {estimated_memory:.2f} GB")

# But the matrix is sparse. Let's check sparsity
# Need to examine the data structure
rows_data = l16_data['rows']
print(f"\nFirst few rows of level-16:")
for i in range(min(3, len(rows_data))):
    row = rows_data[i]
    print(f"  Row {i}: length {len(row)}")
    if isinstance(row, tuple):
        print(f"    First few elements: {row[:4] if len(row) >= 4 else row}")

# Check what format the rows are in
if rows_data:
    first_row = rows_data[0]
    if isinstance(first_row, tuple):
        print(f"\nStructure of first row elements:")
        for j, elem in enumerate(first_row[:4]):
            print(f"  Element {j}: {elem}, type: {type(elem)}")
            if isinstance(elem, tuple) and len(elem) == 2:
                print(f"    word: {elem[0]}, coeff: {elem[1]}, coeff type: {type(elem[1])}")