import sys, importlib, pickle, os, fractions, re, pathlib, sympy as sp
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *

print("="*80)
print("TASK 3: Express vectors in basis and check polynomial dependence on k")
print("="*80)

# Load results
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

degrees = [2, 3, 4, 5]
ranks = {2: 2, 3: 3, 4: 4, 5: 5}

results = {}  # Store results for each degree

for d in degrees:
    print(f"\n--- Degree d = {d} (rank = {ranks[d]}) ---")
    sorted_keys = keys_by_degree[d]
    num_slots = len(poles) * len(sorted_keys)
    
    # Collect vectors and k values
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
    
    # Choose basis: first r = rank[d] vectors
    r = ranks[d]
    basis = vectors[:r]
    basis_k = k_values[:r]
    
    # For each vector (including basis vectors), find coordinates
    # We want to solve: v = sum_{j=0}^{r-1} c_j * basis_j, with c_j in Q(sqrt2)
    
    coordinates = {}  # k -> [c_0, ..., c_{r-1}] where c_j = p_j + q_j*sqrt2
    
    # We'll solve for all vectors at once using linear algebra
    # Convert to sympy for solving
    
    # For each vector v_i, we want c_j^i such that v_i = sum_j c_j^i * basis_j
    # Over Q(sqrt2), c_j^i = p_j^i + q_j^i*sqrt2
    
    # Create matrix equation: V = B * C, where B is basis matrix (over Q(sqrt2))
    # But we need to work over Q by splitting real and sqrt2 parts
    
    # Basis matrix B: each basis vector is length 2*num_slots over Q
    # We want to represent multiplication by (p + q*sqrt2)
    # (p + q*sqrt2)*(a + b*sqrt2) = (p*a + 2*q*b) + (p*b + q*a)*sqrt2
    
    # So for each basis vector j and slot s, with components (a_js, b_js):
    # Contribution to slot s: (p_j*a_js + 2*q_j*b_js, p_j*b_js + q_j*a_js)
    
    # Build big matrix for all vectors
    # We have r unknown pairs (p_j, q_j) for each target vector
    
    # Actually, let's solve for each vector separately
    print(f"  Solving for coordinates...")
    
    all_polynomial = True
    poly_degrees = []
    
    for coord_idx in range(r):  # for each coordinate position
        print(f"\n  Coordinate {coord_idx}:")
        
        # Collect this coordinate value for all k
        coord_values_p = []  # p part (real)
        coord_values_q = []  # q part (sqrt2 part)
        ks_for_fit = []
        
        # We need to solve for each vector to get its coordinates
        # Then check if coord_idx-th coordinate is polynomial in k
        
        # First solve for each vector
        for vec_idx, vec in enumerate(vectors):
            k_val = k_values[vec_idx]
            
            # Solve: vec = sum_{j=0}^{r-1} (p_j + q_j*sqrt2) * basis_j
            # Set up linear system
            
            # Build matrix for this vector
            coeff_rows = []
            target = []
            
            for slot in range(num_slots):
                # Equations for this slot
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
                
                vec_a = vec[2*slot]
                vec_b = vec[2*slot + 1]
                
                coeff_rows.append(row_a)
                coeff_rows.append(row_b)
                target.append(sp.Rational(vec_a.numerator, vec_a.denominator))
                target.append(sp.Rational(vec_b.numerator, vec_b.denominator))
            
            coeff_mat = sp.Matrix(coeff_rows)
            target_vec = sp.Matrix(target)
            
            # Solve linear system
            solution = coeff_mat.solve(target_vec)
            
            # Extract p_{coord_idx} and q_{coord_idx}
            p_val = solution[coord_idx]
            q_val = solution[coord_idx + r]
            
            coord_values_p.append(p_val)
            coord_values_q.append(q_val)
            ks_for_fit.append(k_val)
        
        # Now check if p(k) and q(k) are polynomials in k
        # Use first 5 points to fit polynomial, test on rest
        
        print(f"    Values for k={ks_for_fit[:5]}:")
        print(f"      p = {[str(v) for v in coord_values_p[:5]]}")
        print(f"      q = {[str(v) for v in coord_values_q[:5]]}")
        
        # Try polynomial fits of degree 0..4
        for poly_deg in range(5):
            if len(ks_for_fit) < poly_deg + 1:
                continue
                
            # Fit using first poly_deg+1 points
            fit_ks = ks_for_fit[:poly_deg+1]
            fit_ps = coord_values_p[:poly_deg+1]
            fit_qs = coord_values_q[:poly_deg+1]
            
            # Lagrange interpolation for p
            x = sp.symbols('x')
            
            # Fit polynomial for p
            poly_p = 0
            for i in range(len(fit_ks)):
                term = fit_ps[i]
                for j in range(len(fit_ks)):
                    if i != j:
                        term *= (x - fit_ks[j]) / (fit_ks[i] - fit_ks[j])
                poly_p += term
            
            # Fit polynomial for q
            poly_q = 0
            for i in range(len(fit_ks)):
                term = fit_qs[i]
                for j in range(len(fit_ks)):
                    if i != j:
                        term *= (x - fit_ks[j]) / (fit_ks[i] - fit_ks[j])
                poly_q += term
            
            # Test on all k values
            ok = True
            for idx, k_val in enumerate(ks_for_fit):
                p_exact = coord_values_p[idx]
                q_exact = coord_values_q[idx]
                p_pred = poly_p.subs(x, k_val)
                q_pred = poly_q.subs(x, k_val)
                
                if p_pred != p_exact or q_pred != q_exact:
                    ok = False
                    first_fail = k_val
                    break
            
            if ok:
                print(f"    Found polynomial of degree {poly_deg}:")
                print(f"      p(k) = {sp.simplify(poly_p)}")
                print(f"      q(k) = {sp.simplify(poly_q)}")
                poly_degrees.append(poly_deg)
                break
        
        if len(poly_degrees) <= coord_idx:
            print(f"    No polynomial fit found up to degree 4")
            all_polynomial = False
            poly_degrees.append(None)
    
    print(f"\n  Summary for degree d={d}:")
    if all_polynomial:
        print(f"    All coordinates are polynomials in k")
        print(f"    Polynomial degrees: {poly_degrees}")
    else:
        print(f"    Some coordinates are not polynomials in k")
    
    results[d] = {
        'all_polynomial': all_polynomial,
        'poly_degrees': poly_degrees[:] if all_polynomial else None
    }

print("\n" + "="*80)
print("TASK 3 SUMMARY:")
for d in degrees:
    if results[d]['all_polynomial']:
        print(f"  d={d}: YES, polynomial coordinate functions")
        print(f"       degrees: {results[d]['poly_degrees']}")
    else:
        print(f"  d={d}: NO, coordinate functions not polynomial")

# Save to progress
with open('progress.log', 'a') as f:
    f.write(f"TASK 3 COMPLETE: polynomial check done\n")