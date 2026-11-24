import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import pickle

# Load krows
krows_path = 'research/08_theory/A3f3_level8_compute/k3_newHrows.pickle'
with open(krows_path, 'rb') as f:
    krows_data = pickle.load(f)

print("krows_data type:", type(krows_data))
print("krows_data length:", len(krows_data))
print("\nFirst few entries:")
for i, (k, v) in enumerate(list(krows_data.items())[:5]):
    print(f"  {k}: {v}")