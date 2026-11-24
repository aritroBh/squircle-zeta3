import sys, importlib, pickle, os, fractions, re, pathlib, sympy as sp
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *

def compute_rank_analytical(d):
    """Compute Q(sqrt2)-rank analytically"""
    # Load data for this degree
    all_results = {}
    for k in range(2, 15):
        fn = f"binteg_{k}.pkl"
        with open(fn, 'rb') as f:
            ipf = pickle.load(f)
        res = fastred.reduce_k_fast(k, integrand_pf=ipf)
        all_results[k] = res
    
    AQ = QSI(1,1,0,0)
    BQ = QSI(-1,1,0,0)
    poles = [ZERO, AQ, -AQ, BQ, -BQ]
    
    keys_by_degree = {
        2: [(2, 0), (1, 1)],
        3: [(3, 0), (2, 1)],
        4: [(4, 0), (3, 1), (2, 2)],
        5: [(5, 0), (4, 1), (3, 2)]
    }
    sorted_keys = keys_by_degree[d]
    
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
    
    # Collect vectors for k >= d+1
    vectors = []
    k_values = []
    for k in range(d+1, 15):
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
        
        vectors.append(vec_components)
        k_values.append(k)
    
    # Now compute rank using sympy linear algebra
    # Convert to matrix over Q
    M_q = []
    for vec in vectors:
        M_q.append([sp.Rational(v.numerator, v.denominator) for v in vec])
    
    M = sp.Matrix(M_q)
    rank_q = M.rank()
    
    # Now test if rank over Q(sqrt2) is smaller
    # Try to find if there are sqrt2 relations
    
    # Approach: for candidate rank r, see if we can express
    # all vectors as Q(sqrt2)-linear combinations
    
    for r in range(1, min(len(k_values), 2*num_slots) + 1):
        # Try using first r vectors as basis
        if r > len(vectors):
            break
            
        basis = vectors[:r]
        all_expressible = True
        
        for test_idx in range(r, len(vectors)):
            test_vec = vectors[test_idx]
            
            # Set up equations: test_vec = sum_j (p_j + q_j*sqrt2)*basis[j]
            # Over Q, this gives 2*len(test_vec) equations in 2r unknowns
            
            # We'll solve as linear system over Q
            # Create matrix for the linear system
            
            # For each basis vector j and each slot s:
            # Contribution to slot s: (p_j*a_{j,s} + 2*q_j*b_{j,s}, p_j*b_{j,s} + q_j*a_{j,s})
            
            num_slots = len(poles) * len(sorted_keys)
            num_eqs = 2 * num_slots
            
            # Build coefficient matrix for unknowns p0...p_{r-1}, q0...q_{r-1}
            coeff_matrix = []
            for slot in range(num_slots):
                # Equation for a-component of slot
                row_a = [0]*(2*r)
                # Equation for b-component of slot  
                row_b = [0]*(2*r)
                
                for j in range(r):
                    basis_vec = basis[j]
                    a_js = basis_vec[2*slot]
                    b_js = basis_vec[2*slot + 1]
                    
                    # p_j contributes a_js to a-equation, b_js to b-equation
                    row_a[j] = a_js  # p_j coefficient for a-equation
                    row_b[j + r] = a_js  # q_j coefficient for b-equation (since q_j*a_js)
                    
                    # q_j contributes 2*b_js to a-equation, a_js to b-equation
                    row_a[j + r] = 2 * b_js  # q_j coefficient for a-equation
                    row_b[j] = b_js  # p_j coefficient for b-equation
                
                coeff_matrix.append(row_a)
                coeff_matrix.append(row_b)
            
            # Test vector components
            test_components = []
            for slot in range(num_slots):
                test_components.append(test_vec[2*slot])
                test_components.append(test_vec[2*slot + 1])
            
            # Solve linear system
            coeff_mat = sp.Matrix(coeff_matrix)
            target = sp.Matrix(test_components)
            
            # Check if consistent
            aug = coeff_mat.row_join(target)
            if coeff_mat.rank() != aug.rank():
                all_expressible = False
                break
        
        if all_expressible:
            # Check if basis vectors themselves are Q(sqrt2)-independent
            # Need to check no Q(sqrt2) relation among basis vectors
            basis_independent = True
            
            # Set up: find if there exist (p_j + q_j*sqrt2) not all zero s.t.
            # sum_j (p_j + q_j*sqrt2)*basis_j = 0
            
            # Over Q: for each slot s:
            # sum_j (p_j*a_{j,s} + 2*q_j*b_{j,s}) = 0
            # sum_j (p_j*b_{j,s} + q_j*a_{j,s}) = 0
            
            # This is homogeneous system
            nullspace = coeff_mat.nullspace()
            if nullspace:
                # Check if non-trivial solution exists
                if len(nullspace) > 0:
                    # Actually if nullspace exists, basis vectors are dependent over Q
                    # But we need Q(sqrt2)-dependence
                    # Any Q-nullvector gives Q(sqrt2) relation?
                    # If x = (p0...p_{r-1}, q0...q_{r-1}) is nullvector,
                    # then sum (p_j + q_j*sqrt2)*basis_j = 0
                    basis_independent = False
            
            if basis_independent:
                return r, rank_q
    
    return len(k_values), rank_q  # fallback

print("="*80)
print("TASK 2: Accurate rank computation")
print("="*80)

degrees = [2, 3, 4, 5]
ranks = {}
q_ranks = {}

for d in degrees:
    print(f"\n--- Degree d = {d} ---")
    r, rq = compute_rank_analytical(d)
    ranks[d] = r
    q_ranks[d] = rq
    print(f"  Q(sqrt2)-rank: {r}")
    print(f"  Q-rank: {rq}")

print("\nFinal ranks:")
for d in degrees:
    print(f"  d={d}: Q(sqrt2)-rank = {ranks[d]}, Q-rank = {q_ranks[d]}")
    
    # Check if rank matches expectation
    if ranks[d] == d:
        print(f"    -> Matches d!")
    else:
        print(f"    -> Different from d")