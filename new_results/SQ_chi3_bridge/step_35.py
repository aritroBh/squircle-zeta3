import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import pickle

# Load base rows
base1_path = 'research/08_theory/A3f3_level8_compute/weight3_with_lifts.pickle'
with open(base1_path, 'rb') as f:
    base1_data = pickle.load(f)

print("First row type:", type(base1_data['rows'][0]))
print("First row length:", len(base1_data['rows'][0]))
print("\nFirst row contents:")
for j, (k, v) in enumerate(list(base1_data['rows'][0].items())[:10]):
    print(f"  {k}: {v}")

print("\ncols[:10]:", base1_data['cols'][:10])