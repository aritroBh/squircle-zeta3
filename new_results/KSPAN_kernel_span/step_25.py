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
print("TASK 2: Compute Q(sqrt2)-rank for each degree")
print("="*80)

# First load results
all_results = {}
for k in range(2, 15):
    fn = f"binteg_{k}.pkl"
    with open(fn, 'rb') as f:
        ipf = pickle.load(f)
    res = fastred.reduce_k_fast(k, integrand_pf=ipf)
    all_results[k] = res

print(f"Loaded results for k=2..14")

# Define poles
AQ = QSI(1,1,0,0)  # a = 1+sqrt2
BQ = QSI(-1,1,0,0) # b = sqrt2-1
poles = [ZERO, AQ, -AQ, BQ, -BQ]

def flatten_pf(pf):
    """Convert PF [simple_dict, ...] to vector over Q(sqrt2)"""
    simple_dict = pf[0]
    vec = {}
    for pole in poles:
        if pole in simple_dict:
            coeff = simple_dict[pole]
            vec[pole] = (coeff.p0, coeff.p1)
        else:
            vec[pole] = (fractions.Fraction(0), fractions.Fraction(0))
    return vec

degrees = [2, 3, 4, 5]
keys_by_degree = {
    2: [(2, 0), (1, 1)],
    3: [(3, 0), (2, 1)],
    4: [(4, 0), (3, 1), (2, 2)],
    5: [(5, 0), (4, 1), (3, 2)]
}

# Compute rank for each degree using sympy DomainMatrix
ranks = {}
bases_vectors = {}  # Store basis vectors for each degree
bases_indices = {}  # Store which k values form a basis

for d in degrees:
    print(f"\n--- Degree d = {d} ---")
    sorted_keys = keys_by_degree[d]
    
    # Collect vectors for k >= d+1
    vectors = []
    k_values = []
    for k in range(d+1, 15):
        if k not in all_results:
            continue
        
        vec_components = []
        for key in sorted_keys:
            if key in all_results[k]:
                pf = all_results[k][key]
                vec_key = flatten_pf(pf)
                for pole in poles:
                    p0, p1 = vec_key[pole]
                    vec_components.append(sp.Rational(p0.numerator, p0.denominator))
                    vec_components.append(sp.Rational(p1.numerator, p1.denominator))
            else:
                for _ in range(len(poles) * 2):
                    vec_components.append(sp.Rational(0))
        
        vectors.append(vec_components)
        k_values.append(k)
    
    # Create matrix: rows are vectors for different k, columns are components
    # Each component is in Q(sqrt2), represented as a + b*sqrt2
    # We'll create a matrix over Q(sqrt2) algebraic field
    
    # First, let's compute rank by considering each component as independent over Q
    # Actually, we can treat (a + b*sqrt2) as vector (a, b) over Q
    # So we'll convert each Q(sqrt2) component to two rational components
    
    num_components = len(poles) * 2 * len(sorted_keys)  # each pole has 2 components
    # Actually we already doubled: each Q(sqrt2) component is already split into p0 and p1
    # So num_components is already the total
    
    # Convert to sympy matrix over Q
    mat_rows = []
    for vec in vectors:
        # vec is already list of Fractions, convert to Rationals
        mat_rows.append([sp.Rational(v.numerator, v.denominator) for v in vec])
    
    M = sp.Matrix(mat_rows)
    print(f"  Matrix shape: {M.shape} (rows={len(k_values)}, cols={num_components})")
    
    # Compute rank over Q
    rank_q = M.rank()
    
    # Now we need to interpret this as Q(sqrt2) rank
    # If vectors are over Q(sqrt2), each Q(sqrt2) dimension corresponds to 2 Q dimensions
    # But careful: our splitting already separated real and sqrt2 parts
    # Actually, if we have r Q(sqrt2)-independent vectors, that means 2r Q-independent vectors
    # But the converse isn't necessarily true...
    
    # Better: create algebraic field Q(sqrt2) and construct matrix directly
    sqrt2 = sp.sqrt(2)
    Qsqrt2 = sp.QQ.algebraic_field(sqrt2)
    
    # Convert vectors to Q(sqrt2) elements
    # Each pair (p0, p1) represents p0 + p1*sqrt2
    # But we already separated them... Let me reconstruct
    
    # Actually, for each key and pole, we have two components in our vector
    # We want to combine them back into Q(sqrt2) elements
    
    # Reconstruct: for each key-pole pair, combine p0 and p1
    mat_rows_qsqrt2 = []
    for vec in vectors:
        row = []
        # vec has structure: for each key, for each pole: (p0, p1)
        # So we have groups of 2
        for i in range(0, len(vec), 2):
            p0 = vec[i]
            p1 = vec[i+1]
            # Create element in Q(sqrt2)
            elem = sp.Rational(p0.numerator, p0.denominator) + sp.Rational(p1.numerator, p1.denominator) * sqrt2
            row.append(elem)
        mat_rows_qsqrt2.append(row)
    
    # Create domain matrix over Q(sqrt2)
    try:
        # Try using DomainMatrix
        from sympy.polys.matrices import DomainMatrix
        
        # Convert to list of lists
        rows_data = [[elem for elem in row] for row in mat_rows_qsqrt2]
        
        # Create DomainMatrix
        dm = DomainMatrix.from_list(rows_data, (len(k_values), len(mat_rows_qsqrt2[0])), Qsqrt2)
        
        # Compute rank
        rank_qsqrt2 = dm.rank()
        print(f"  Rank over Q(sqrt2): {rank_qsqrt2}")
        
        # Also compute rank over Q by thinking about it
        print(f"  Rank over Q (from matrix): {rank_q}")
        print(f"  Note: over Q, max possible rank would be {2*rank_qsqrt2}")
        
        # Store results
        ranks[d] = rank_qsqrt2
        
        # Find basis: first 'rank_qsqrt2' linearly independent vectors
        # We'll use sympy to find pivot columns
        M_qsqrt2 = sp.Matrix(mat_rows_qsqrt2)
        # Get row echelon form
        M_rref, pivots = M_qsqrt2.rref()
        print(f"  Pivot columns (0-indexed): {pivots}")
        
        # Basis vectors are those corresponding to pivot columns
        basis_k_indices = []
        if len(pivots) > 0:
            # The pivots tell us which original vectors (rows) correspond to independent directions
            # Actually pivots are column indices, we need to think differently
            # For basis, we need a set of k values whose vectors are independent
            
            # Approach: start with empty set, add vectors greedily
            basis_vectors = []
            basis_k = []
            independent_set = []
            
            for idx, k_val in enumerate(k_values):
                vec = mat_rows_qsqrt2[idx]
                # Represent as DomainMatrix row
                temp_vectors = basis_vectors + [vec]
                # Check if linear independent
                temp_dm = DomainMatrix.from_list(temp_vectors, (len(temp_vectors), len(vec)), Qsqrt2)
                if temp_dm.rank() == len(temp_vectors):
                    basis_vectors.append(vec)
                    basis_k.append(k_val)
                    if len(basis_vectors) >= rank_qsqrt2:
                        break
            
            print(f"  Basis k values: {basis_k}")
            bases_indices[d] = basis_k
            bases_vectors[d] = basis_vectors
        
    except Exception as e:
        print(f"  Error computing rank with DomainMatrix: {e}")
        print(f"  Falling back to Q-rank analysis")
        
        # Alternative: each Q(sqrt2) component gives 2 Q dimensions
        # If rank over Q is r, then Q(sqrt2) rank is at least ceil(r/2)
        # But could be less if there are sqrt2 relations
        
        # Let's create matrix over Q with each Q(sqrt2) component as 2 columns
        # Actually we already have M which is over Q with separated components
        
        # Try to see if pairs (p0, p1) for each component are related
        # For each "slot" (key-pole position), we have pairs (p0, p1) across different k
        # Create matrix for each slot and check rank
        
        slot_rank = 0
        num_slots = len(poles) * len(sorted_keys)
        print(f"  Number of Q(sqrt2) slots: {num_slots}")
        
        for slot_idx in range(num_slots):
            # Extract p0 and p1 for this slot across all k
            p0_col = []
            p1_col = []
            for vec in vectors:
                p0 = vec[2*slot_idx]
                p1 = vec[2*slot_idx + 1]
                p0_col.append(sp.Rational(p0.numerator, p0.denominator))
                p1_col.append(sp.Rational(p1.numerator, p1.denominator))
            
            # Check if they're linearly independent over Q
            # Create matrix of these two columns
            slot_mat = sp.Matrix([p0_col, p1_col]).T  # rows are k values, columns are p0, p1
            slot_rank += slot_mat.rank()
        
        print(f"  Sum of slot ranks over Q: {slot_rank}")
        # Actually this doesn't give Q(sqrt2) rank directly...
        
        # Let me try a simpler heuristic: use first few k values
        ranks[d] = "TODO"
        continue
    
print("\nSummary of ranks:")
for d in degrees:
    if d in ranks:
        print(f"  d={d}: rank = {ranks[d]}")

# Progress log
with open('progress.log', 'a') as f:
    f.write(f"TASK 2 in progress: rank computation attempted\n")