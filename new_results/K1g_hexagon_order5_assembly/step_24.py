import os, pickle
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f, 'E5_polys.pkl'), 'rb') as f:
    E5_polys = pickle.load(f)
E5 = E5_polys[5]

# Separate local (c=0) and non-local (c>=1)
local_keys = [k for k in E5 if k[3]==0]
nonlocal_keys = [k for k in E5 if k[3]>=1]
print(f"E5 total keys: {len(E5)}")
print(f"local (c=0): {len(local_keys)}")
print(f"non-local (c>=1): {len(nonlocal_keys)}")
print()
# Check what symbols appear in coefficients
import sympy as sp
all_syms = set()
for k,v in E5.items():
    all_syms |= v.free_symbols
print("All free symbols in E5 coefficients:")
print(sorted([str(s) for s in all_syms]))
