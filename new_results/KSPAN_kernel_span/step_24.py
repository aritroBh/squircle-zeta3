import sys, importlib, time, pickle, os, json, fractions, re, pathlib
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *
import buildintegrand
importlib.reload(buildintegrand)
from buildintegrand import build_integrand_pf
import sympy as sp

print("="*80)
print("TASK 2: Flatten degree-d pieces into vectors and compute ranks")
print("="*80)

# First load all results from task 1 (or recompute)
all_results = {}
for k in range(2, 15):
    fn = f"binteg_{k}.pkl"
    with open(fn, 'rb') as f:
        ipf = pickle.load(f)
    res = fastred.reduce_k_fast(k, integrand_pf=ipf)
    all_results[k] = res

print(f"Loaded results for k=2..14")

# Define poles of interest: 0, a, -a, b, -b where a=1+sqrt2, b=sqrt2-1
AQ = QSI(1,1,0,0)  # a = 1+sqrt2
BQ = QSI(-1,1,0,0) # b = sqrt2-1
poles = [ZERO, AQ, -AQ, BQ, -BQ]

# Helper to flatten a PF into vector
def flatten_pf(pf):
    """Convert PF [simple_dict, ...] to vector over Q(sqrt2)"""
    simple_dict = pf[0]
    # Create vector indexed by pole order: for each pole get (p0, p1)
    # We'll represent as tuple (p0_rational, p1_rational)
    vec = {}
    for pole in poles:
        if pole in simple_dict:
            coeff = simple_dict[pole]
            vec[pole] = (coeff.p0, coeff.p1)  # (real part, sqrt2 part)
        else:
            vec[pole] = (fractions.Fraction(0), fractions.Fraction(0))
    return vec

# Now for each degree d = 2,3,4,5, collect vectors for k >= d+1
degrees = [2, 3, 4, 5]
vectors_by_degree = {}
keys_by_degree = {}

for d in degrees:
    print(f"\n--- Degree d = {d} ---")
    
    # Find all (i,j) keys with i+j = d
    all_keys = set()
    for k in all_results:
        for key in all_results[k]:
            i, j = key
            if i + j == d:
                all_keys.add(key)
    
    # Sort keys: (i,j) with i >= j, weight = i+j+1
    sorted_keys = sorted(all_keys, key=lambda x: (-x[0], x[1]))  # sort by i descending, then j
    keys_by_degree[d] = sorted_keys
    print(f"  Keys with i+j={d}: {sorted_keys}")
    
    # For each k >= d+1, flatten all degree-d pieces into one big vector
    vectors = {}
    for k in range(d+1, 15):  # k >= d+1
        if k not in all_results:
            continue
            
        vec_components = []
        for key in sorted_keys:
            if key in all_results[k]:
                pf = all_results[k][key]
                vec_key = flatten_pf(pf)
                # Flatten further: for each pole, add p0 and p1 to vector
                for pole in poles:
                    p0, p1 = vec_key[pole]
                    vec_components.append(p0)
                    vec_components.append(p1)
            else:
                # Key missing, add zeros
                for _ in range(len(poles) * 2):  # 2 components per pole
                    vec_components.append(fractions.Fraction(0))
        
        vectors[k] = vec_components
    
    vectors_by_degree[d] = vectors
    
    print(f"  Vectors collected for k={list(vectors.keys())}")
    if vectors:
        vec_len = len(next(iter(vectors.values())))
        print(f"  Vector length: {vec_len}")
        
        # Let's examine a small sample
        k_sample = d+1
        if k_sample in vectors:
            sample_vec = vectors[k_sample]
            print(f"  Sample vector for k={k_sample} (first 10 components): {sample_vec[:10]}")

print("\nDone collecting vectors")