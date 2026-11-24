import sys, importlib, pickle, os, re, pathlib, fractions
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *
import sympy as sp

S2 = sp.sqrt(2)

# Let me examine k=3 (3,0) more carefully
k = 3
fn = f"binteg_{k}.pkl"
with open(fn, 'rb') as f:
    ipf = pickle.load(f)

res = fastred.reduce_k_fast(k, integrand_pf=ipf)
pf = res[(3,0)]
print(f"k={k}, key=(3,0):")
print(f"  PF simple dict: {pf[0]}")

# Parse the reference file directly
outputs_dir = pathlib.Path("<project>/paper/squircle_note/outputs")
k123_file = outputs_dir / "ibp_reduce_k123.out"
with open(k123_file, 'r') as f:
    content = f.read()

# Find k=3 section
lines = content.split('\n')
in_k3 = False
for i, line in enumerate(lines):
    if line.strip() == "===== k = 3:":
        in_k3 = True
    elif line.strip().startswith("===== k ="):
        if in_k3:
            break
    if in_k3 and "l1^3 l2^0" in line:
        print(f"\nReference lines for (3,0):")
        j = i + 1
        while j < len(lines) and lines[j].strip() and not lines[j].strip().startswith("weight"):
            print(f"  {lines[j].strip()}")
            j += 1
        break