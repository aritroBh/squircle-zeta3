import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
from fractions import Fraction
import pickle

# Load base rows
base1_path = 'research/08_theory/A3f3_level8_compute/weight3_with_lifts.pickle'
with open(base1_path, 'rb') as f:
    base1_data = pickle.load(f)

print("base1_data keys:", base1_data.keys())
print("cols type:", type(base1_data['cols']))
print("rows type:", type(base1_data['rows']))
print("cols length:", len(base1_data['cols']))
print("rows length:", len(base1_data['rows']))

# Check what's in rows
print("\nFirst few rows:")
for i, row in enumerate(list(base1_data['rows'].values())[:3]):
    print(f"Row {i}: {len(row)} entries")
    # Show some entries
    for j, (k, v) in enumerate(list(row.items())[:5]):
        print(f"  {k}: {v}")