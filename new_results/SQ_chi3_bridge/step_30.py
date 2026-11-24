import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import pickle

# Load vals
vals_path = 'research/08_theory/A3f6_level16/vals8_150digit.pickle'
with open(vals_path, 'rb') as f:
    vals = pickle.load(f)

# Check the structure - what words are in vals?
print("First 20 words in vals:")
for i, w in enumerate(sorted(vals.keys())[:20]):
    print(f"  {w}")

print(f"\nTotal words in vals: {len(vals)}")

# Check what Lv contains
Lv = C.atom(('log', C.A))
print(f"\nLv words: {list(Lv.keys())}")