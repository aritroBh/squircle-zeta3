import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import pickle

# Load base rows
base1_path = 'research/08_theory/A3f3_level8_compute/weight3_with_lifts.pickle'
with open(base1_path, 'rb') as f:
    base1_data = pickle.load(f)

print("First row:")
print(base1_data['rows'][0])

print("\nSecond row:")
print(base1_data['rows'][1])

print("\nrank_Q:", base1_data['rank_Q'])