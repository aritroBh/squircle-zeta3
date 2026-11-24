import sys, importlib, pickle, os, fractions, re, pathlib, sympy as sp
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *

print("="*80)
print("TASK 2: Compute rank using rational representation")
print("="*80)

# Load results
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

# Compute rank by treating (a + b*sqrt2) as 2D vector (a,b) over Q
# The Q(sqrt2)-rank r corresponds to rank over Q being 2r, with structure
ranks = {}

for d in degrees:
    print(f"\n--- Degree d = {d} ---")
    sorted_keys = keys_by_degree[d]
    
    # Collect vectors for k >= d+1
    mat_rows = []
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
                    vec_components.append(p0)
                    vec_components.append(p1)
            else:
                for _ in range(len(poles) * 2):
                    vec_components.append(fractions.Fraction(0))
        
        mat_rows.append([sp.Rational(v.numerator, v.denominator) for v in vec_components])
        k_values.append(k)
    
    M = sp.Matrix(mat_rows)
    print(f"  Matrix shape: {M.shape} (rows={len(k_values)}, cols={len(mat_rows[0])})")
    
    # Compute rank over Q
    rank_q = M.rank()
    print(f"  Rank over Q: {rank_q}")
    
    # Now we need Q(sqrt2) rank
    # Each Q(sqrt2) dimension contributes 2 to rank over Q if the sqrt2 part is independent
    # But there could be linear dependencies involving sqrt2
    
    # Create matrix for each "slot" (key-pole position)
    # Each slot has structure: across k, we have pairs (a_k, b_k)
    # representing a_k + b_k*sqrt2
    
    num_slots = len(poles) * len(sorted_keys)
    slot_data = []
    
    for slot_idx in range(num_slots):
        a_col = []
        b_col = []
        for row in mat_rows:
            a_col.append(row[2*slot_idx])
            b_col.append(row[2*slot_idx + 1])
        slot_data.append((a_col, b_col))
    
    # Now analyze possible sqrt2 relations
    # A Q(sqrt2) relation would be: sum_i c_i*(a_i + b_i*sqrt2) = 0, c_i in Q(sqrt2)
    # Write c_i = p_i + q_i*sqrt2
    # Then: sum_i p_i*a_i + 2*q_i*b_i = 0 and sum_i p_i*b_i + q_i*a_i = 0
    
    # Instead, let's try to find Q(sqrt2)-rank directly
    # Approach: look for linear dependencies over Q(sqrt2) by solving
    # M * x = 0 where x has elements in Q(sqrt2)
    
    # Convert to complex representation: treat sqrt2 as algebraic number
    # Actually symmetric: map a + b*sqrt2 to vector (a, b) but with multiplication rules
    
    # Simpler: Since we suspect small rank, let's check if first few k values are independent
    # and if others are in their span
    
    # Try to express vectors as linear combinations with Q(sqrt2) coefficients
    # We'll test small candidate ranks
    
    max_possible_rank = min(len(k_values), num_slots)
    
    for r in range(1, max_possible_rank + 1):
        # Try using first r vectors as basis
        basis = mat_rows[:r]
        
        # Try to express all vectors in this basis with Q(sqrt2) coefficients
        all_in_span = True
        
        for test_idx, test_vec in enumerate(mat_rows):
            # Skip basis vectors themselves
            if test_idx < r:
                continue
                
            # Try to solve: test_vec = sum_{j=0}^{r-1} (p_j + q_j*sqrt2) * basis[j]
            # Over Q, this is 2r unknowns (p_j, q_j) and len(test_vec) equations
            
            # Set up linear system over Q
            # For each component i: test_vec[i] = sum_j (p_j * basis[j][i] + 2*q_j * basis[j][i] for sqrt2 part?)
            # Actually careful: (p + q*sqrt2)*(a + b*sqrt2) = (p*a + 2*q*b) + (p*b + q*a)*sqrt2
            
            # So for each slot (pair of components), we have:
            # test_a = sum_j (p_j*a_j + 2*q_j*b_j)
            # test_b = sum_j (p_j*b_j + q_j*a_j)
            
            # We have 2*num_slots equations, 2r unknowns
            
            # Let's implement this
            unknowns = []
            equations = []
            
            # Create variables p0...p_{r-1}, q0...q_{r-1}
            p_vars = [sp.symbols(f'p{j}') for j in range(r)]
            q_vars = [sp.symbols(f'q{j}') for j in range(r)]
            unknowns = p_vars + q_vars
            
            # For each slot
            for slot_idx in range(num_slots):
                test_a = test_vec[2*slot_idx]
                test_b = test_vec[2*slot_idx + 1]
                
                expr_a = 0
                expr_b = 0
                for j in range(r):
                    basis_a = basis[j][2*slot_idx]
                    basis_b = basis[j][2*slot_idx + 1]
                    expr_a += p_vars[j] * basis_a + 2 * q_vars[j] * basis_b
                    expr_b += p_vars[j] * basis_b + q_vars[j] * basis_a
                
                equations.append(sp.Eq(expr_a, test_a))
                equations.append(sp.Eq(expr_b, test_b))
            
            # Solve linear system
            try:
                solution = sp.solve(equations, unknowns, dict=True)
                if not solution:
                    all_in_span = False
                    break
            except:
                all_in_span = False
                break
        
        if all_in_span:
            print(f"  Candidate Q(sqrt2)-rank r = {r}: all vectors in span of first {r}")
            ranks[d] = r
            break
    
    if d not in ranks:
        print(f"  Could not determine Q(sqrt2)-rank, using min({len(k_values)}, {num_slots})")
        ranks[d] = min(len(k_values), num_slots)

print("\nComputed ranks:")
for d in degrees:
    print(f"  d={d}: rank = {ranks[d]}")

# Save to progress log
with open('progress.log', 'a') as f:
    f.write(f"TASK 2: Computed ranks: {ranks}\n")