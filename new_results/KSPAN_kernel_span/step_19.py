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

S2 = sp.sqrt(2)

def parse_rational_expr(expr_str):
    """Parse expression like 'sqrt(2)/(16*v)' or '(2*sqrt(2) + 7)/(16*v)'"""
    expr_str = expr_str.strip()
    
    try:
        expr = sp.sympify(expr_str)
        expr = sp.expand(expr)
        
        # Get numerator and denominator
        num, den = expr.as_numer_denom()
        
        # Check denominator format
        if den == sp.Symbol('v'):
            pole = ZERO
            coeff = num
        elif str(den).startswith('('):
            # Pattern: (v +/- something)
            den_str = str(den)
            if den_str.startswith('(') and den_str.endswith(')'):
                den_str = den_str[1:-1]
            
            # Solve v - pole = 0
            v = sp.Symbol('v')
            pole_expr = sp.solve(sp.Eq(den, 0), v)[0]
            
            # Should be real
            pole_real = sp.re(pole_expr)
            pole_imag = sp.im(pole_expr)
            assert pole_imag == 0
            
            # Convert pole_real to Q(sqrt2)
            p0 = sp.simplify(pole_real.subs(S2, 0))
            p1 = sp.simplify((pole_real - p0)/S2)
            
            pole = QSI(fractions.Fraction(int(sp.Rational(p0).p), int(sp.Rational(p0).q)),
                      fractions.Fraction(int(sp.Rational(p1).p), int(sp.Rational(p1).q)))
            coeff = num
        else:
            # Try other patterns
            if str(den).endswith('*v'):
                # Pattern like 16*v
                coeff_const = sp.simplify(sp.S(str(den).replace('*v', '')))
                pole = ZERO
                coeff = sp.simplify(num / coeff_const)
            else:
                # Maybe it's something like "8*v + 8 + 8*sqrt(2)"
                # Try solving directly
                v = sp.Symbol('v')
                pole_expr = sp.solve(sp.Eq(den, 0), v)[0]
                pole_real = sp.re(pole_expr)
                pole_imag = sp.im(pole_expr)
                assert pole_imag == 0
                
                p0 = sp.simplify(pole_real.subs(S2, 0))
                p1 = sp.simplify((pole_real - p0)/S2)
                
                pole = QSI(fractions.Fraction(int(sp.Rational(p0).p), int(sp.Rational(p0).q)),
                          fractions.Fraction(int(sp.Rational(p1).p), int(sp.Rational(p1).q)))
                coeff = num
        
        # Convert coeff to QSI
        coeff = sp.expand(coeff)
        c0 = sp.simplify(coeff.subs(S2, 0))
        c1 = sp.simplify((coeff - c0)/S2)
        
        coeff_qsi = QSI(fractions.Fraction(int(sp.Rational(c0).p), int(sp.Rational(c0).q)),
                       fractions.Fraction(int(sp.Rational(c1).p), int(sp.Rational(c1).q)))
        
        return pole, coeff_qsi
        
    except Exception as e:
        print(f"Error parsing '{expr_str}': {e}")
        raise

# Define Rw2
AQ = QSI(1,1,0,0)  # a = 1+sqrt2
BQ = QSI(-1,1,0,0) # b = sqrt2-1
Rw2pf = [{ZERO: QSI(fractions.Fraction(1,2)), 
          AQ: QSI(fractions.Fraction(1,2)), 
          -AQ: QSI(-1),
          -BQ: QSI(fractions.Fraction(-1,2))},
         [], [], ZERO]

def parse_ref_section(content, k_value):
    """Parse reference output section for given k"""
    sections = {}
    
    # Find sections like "===== k = 3:" 
    lines = content.split('\n')
    current_k = None
    current_weight = None
    current_key = None
    current_pf = {}
    
    for line in lines:
        line = line.strip()
        
        # Match k section header
        k_match = re.match(r'=====\s*k\s*=\s*(\d+):', line)
        if k_match:
            current_k = int(k_match.group(1))
            if current_k == k_value:
                # Reset for new k
                sections = {}
            continue
        
        if current_k != k_value:
            continue
            
        # Match weight line
        weight_match = re.match(r'weight\s+(\d+)\s+l1\^(\d+)\s+l2\^(\d+):', line)
        if weight_match:
            weight = int(weight_match.group(1))
            i = int(weight_match.group(2))
            j = int(weight_match.group(3))
            current_key = (i, j)
            current_weight = weight
            current_pf = {}
            sections[current_key] = current_pf
            continue
        
        # Parse rational expressions for current key
        if current_key and line and not line.startswith('====='):
            # Skip empty lines
            if re.match(r'^\s*$', line):
                continue
                
            # Parse the rational expression
            try:
                pole, coeff = parse_rational_expr(line)
                if pole in current_pf:
                    current_pf[pole] = current_pf[pole] + coeff
                else:
                    current_pf[pole] = coeff
            except Exception as e:
                print(f"Warning: Could not parse '{line}': {e}")
    
    return sections

def compare_pf(computed_pf, ref_pf_dict):
    """Compare computed PF [simple_dict, ...] with reference dict {pole: coeff}"""
    # computed_pf is [simple_dict, higher_plus_i, higher_minus_i, const]
    computed_simple = computed_pf[0]
    
    # Check that computed has no higher parts or const
    if computed_pf[1] or computed_pf[2] or not computed_pf[3].is_zero():
        return False, "Has higher parts or non-zero const"
    
    # Compare simple poles
    if len(computed_simple) != len(ref_pf_dict):
        return False, f"Different number of poles: {len(computed_simple)} vs {len(ref_pf_dict)}"
    
    # Compare each pole
    for pole, coeff in ref_pf_dict.items():
        if pole not in computed_simple:
            return False, f"Missing pole {pole}"
        if computed_simple[pole] != coeff:
            return False, f"Different coeff at pole {pole}: {computed_simple[pole]} vs {coeff}"
    
    return True, "OK"

# Now run task 1
print("="*80)
print("TASK 1: Running reduce_k_fast for k=2..14 and checking controls")
print("="*80)

outputs_dir = pathlib.Path("<project>/paper/squircle_note/outputs")
k123_file = outputs_dir / "ibp_reduce_k123.out"
k4_file = outputs_dir / "ibp_reduce_k4.out"

# Parse reference files
k3_ref = {}
k4_ref = {}

if k123_file.exists():
    with open(k123_file, 'r') as f:
        k123_content = f.read()
    k3_ref = parse_ref_section(k123_content, 3)

if k4_file.exists():
    with open(k4_file, 'r') as f:
        k4_content = f.read()
    k4_ref = parse_ref_section(k4_content, 4)

print(f"\nParsed reference keys:")
print(f"  k=3: {list(k3_ref.keys())}")
print(f"  k=4: {list(k4_ref.keys())}")
if k3_ref:
    print(f"\nSample from k=3 (2,1): {list(k3_ref[(2,1)].items())[:3]}")
if k4_ref:
    print(f"Sample from k=4 (2,1): {list(k4_ref[(2,1)].items())[:3]}")

# Now run reduction for k=2..14
all_results = {}
all_controls_passed = True

for k in range(2, 15):
    print(f"\n--- k = {k} ---")
    fn = f"binteg_{k}.pkl"
    
    with open(fn, 'rb') as f:
        ipf = pickle.load(f)
    
    res = fastred.reduce_k_fast(k, integrand_pf=ipf)
    all_results[k] = res
    
    # Control 1: (1,0) == Rw2
    d = pf_sub(res[(1,0)], Rw2pf)
    kernel_eq = (not d[0]) and (not d[1]) and (not d[2]) and d[3].is_zero()
    print(f"  (1,0) == Rw2: {kernel_eq}")
    if not kernel_eq:
        all_controls_passed = False
    
    # Control 2: All PFs have no higher parts, zero const, coefficients in Q(sqrt2)
    all_good = True
    for key, pf in res.items():
        if pf[1] or pf[2]:
            all_good = False
            print(f"  ERROR: PF {key} has higher parts")
            break
        if not pf[3].is_zero():
            all_good = False
            print(f"  ERROR: PF {key} has non-zero const")
            break
        for pole, coeff in pf[0].items():
            if coeff.p2 != 0 or coeff.p3 != 0:
                all_good = False
                print(f"  ERROR: PF {key} pole {pole} has complex coeff: {coeff}")
                break
        if not all_good:
            break
    
    print(f"  All PFs good: {all_good}")
    if not all_good:
        all_controls_passed = False
    
    # Control 3: Compare with reference for k=3 and k=4 degree-2 keys
    if k == 3:
        print(f"\n  Comparing degree-2 keys with reference k=3:")
        for key in [(2,1), (3,0)]:  # degree-2 keys for k=3
            if key in res and key in k3_ref:
                ok, msg = compare_pf(res[key], k3_ref[key])
                print(f"    {key}: {ok} ({msg})")
                if not ok:
                    all_controls_passed = False
            else:
                print(f"    {key}: Missing in computed or reference")
                all_controls_passed = False
    
    if k == 4:
        print(f"\n  Comparing degree-2 keys with reference k=4:")
        for key in [(2,1), (3,0)]:  # degree-2 keys for k=4
            if key in res and key in k4_ref:
                ok, msg = compare_pf(res[key], k4_ref[key])
                print(f"    {key}: {ok} ({msg})")
                if not ok:
                    all_controls_passed = False
            else:
                print(f"    {key}: Missing in computed or reference")
                all_controls_passed = False

print("\n" + "="*80)
print(f"TASK 1 SUMMARY: All controls passed: {all_controls_passed}")
print("="*80)

if all_controls_passed:
    # Save progress
    with open('progress.log', 'a') as f:
        f.write(f"TASK 1 COMPLETE: All controls passed for k=2..14\n")
else:
    raise ValueError("Some controls failed")