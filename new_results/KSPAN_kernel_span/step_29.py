import sys, importlib, pickle, os, fractions, re, pathlib, sympy as sp
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *

print("="*80)
print("TASK 2: Compute exact Q(sqrt2)-rank")
print("="*80)

# Load all results once
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

keys_by_degree = {
    2: [(2, 0), (1, 1)],
    3: [(3, 0), (2, 1)],
    4: [(4, 0), (3, 1), (2, 2)],
    5: [(5, 0), (4, 1), (3, 2)]
}

def compute_degree_rank(d):
    sorted_keys = keys_by_degree[d]
    num_slots = len(poles) * len(sorted_keys)
    
    # Collect vectors
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
    
    # Convert to sympy Rationals
    M_q = []
    for vec in vectors:
        M_q.append([sp.Rational(v.numerator, v.denominator) for v in vec])
    
    M = sp.Matrix(M_q)
    rank_q = M.rank()
    
    # Now find Q(sqrt2)-rank
    # We'll test if vectors are in Q(sqrt2)-span of first r vectors
    max_r = min(len(k_values), num_slots)
    
    for r in range(1, max_r + 1):
        if r > len(vectors):
            break
            
        # Check if first r vectors are Q(sqrt2)-independent
        # Build matrix for basis independence check
        basis = vectors[:r]
        
        # Check independence: are there p_j, q_j in Q (not all zero) s.t.
        # sum_j (p_j + q_j*sqrt2)*basis_j = 0?
        # Over Q: for each slot s:
        # sum_j (p_j*a_{j,s} + 2*q_j*b_{j,s}) = 0
        # sum_j (p_j*b_{j,s} + q_j*a_{j,s}) = 0
        
        # Build coefficient matrix for homogeneous system
        coeff_rows = []
        for slot in range(num_slots):
            row_a = [0]*(2*r)
            row_b = [0]*(2*r)
            
            for j in range(r):
                basis_vec = basis[j]
                a_js = basis_vec[2*slot]
                b_js = basis_vec[2*slot + 1]
                
                # p_j contributes
                row_a[j] = sp.Rational(a_js.numerator, a_js.denominator)
                row_b[j + r] = sp.Rational(a_js.numerator, a_js.denominator)
                
                # q_j contributes
                row_a[j + r] = sp.Rational(2*b_js.numerator, b_js.denominator)
                row_b[j] = sp.Rational(b_js.numerator, b_js.denominator)
            
            coeff_rows.append(row_a)
            coeff_rows.append(row_b)
        
        coeff_mat = sp.Matrix(coeff_rows)
        nullspace_dim = coeff_mat.shape[1] - coeff_mat.rank()
        
        if nullspace_dim > 0:
            # Has non-trivial solutions
            continue  # basis dependent, try larger r
        
        # Now check if all other vectors are in Q(sqrt2)-span
        all_in_span = True
        for test_idx in range(r, len(vectors)):
            test_vec = vectors[test_idx]
            
            # Build augmented system for this test vector
            aug_rows = []
            for slot in range(num_slots):
                row_a = [0]*(2*r)
                row_b = [0]*(2*r)
                
                for j in range(r):
                    basis_vec = basis[j]
                    a_js = basis_vec[2*slot]
                    b_js = basis_vec[2*slot + 1]
                    
                    row_a[j] = sp.Rational(a_js.numerator, a_js.denominator)
                    row_a[j + r] = sp.Rational(2*b_js.numerator, b_js.denominator)
                    
                    row_b[j] = sp.Rational(b_js.numerator, b_js.denominator)
                    row_b[j + r] = sp.Rational(a_js.numerator, a_js.denominator)
                
                test_a = test_vec[2*slot]
                test_b = test_vec[2*slot + 1]
                
                aug_rows.append(row_a + [sp.Rational(test_a.numerator, test_a.denominator)])
                aug_rows.append(row_b + [sp.Rational(test_b.numerator, test_b.denominator)])
            
            aug_mat = sp.Matrix(aug_rows)
            coeff_part = aug_mat[:, :-1]
            target_part = aug_mat[:, -1]
            
            # Check consistency
            if coeff_part.rank() != aug_mat.rank():
                all_in_span = False
                break
        
        if all_in_span:
            return r, rank_q
    
    return max_r, rank_q

degrees = [2, 3, 4, 5]
ranks = {}
q_ranks = {}

for d in degrees:
    print(f"\n--- Degree d = {d} ---")
    r, rq = compute_degree_rank(d)
    ranks[d] = r
    q_ranks[d] = rq
    print(f"  Q(sqrt2)-rank: {r}")
    print(f"  Q-rank: {rq}")

print("\n" + "="*80)
print("RESULTS:")
for d in degrees:
    print(f"  d={d}: Q(sqrt2)-rank = {ranks[d]}, Q-rank = {q_ranks[d]}")

# Save to progress
with open('progress.log', 'a') as f:
    f.write(f"TASK 2 COMPLETE: ranks = {ranks}\n")