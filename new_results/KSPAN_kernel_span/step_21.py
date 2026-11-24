import sys, importlib, re, sympy as sp, fractions
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *

S2 = sp.sqrt(2)

# Test parsing the specific expression
expr_str = "-1/(3*v + 3 + 3*sqrt(2))"
print(f"Parsing: {expr_str}")

try:
    expr = sp.sympify(expr_str)
    print(f"  sympified: {expr}")
    expr = sp.expand(expr)
    print(f"  expanded: {expr}")
    
    num, den = expr.as_numer_denom()
    print(f"  numerator: {num}")
    print(f"  denominator: {den}")
    
    # Solve for pole
    v = sp.Symbol('v')
    pole_expr = sp.solve(sp.Eq(den, 0), v)[0]
    print(f"  pole expr: {pole_expr}")
    
    # Should be: -1 - sqrt(2)
    pole_real = sp.re(pole_expr)
    pole_imag = sp.im(pole_expr)
    print(f"  pole real: {pole_real}")
    
    # Convert to Q(sqrt2)
    p0 = sp.simplify(pole_real.subs(S2, 0))
    p1 = sp.simplify((pole_real - p0)/S2)
    print(f"  p0: {p0}, p1: {p1}")
    
    # The coefficient should be num * (derivative of denominator)
    # Actually for 1/(v - pole), coefficient is 1
    # But here we have 1/[3*(v - pole)], coefficient is 1/3
    # My parser is just extracting num = -1
    # I need to account for the factor in denominator!
    
    # Denominator is 3*(v - pole), so coefficient should be -1/3
    # I need to find the constant factor in denominator
    
    # Factor out the constant
    den_factored = sp.factor(den)
    print(f"  denominator factored: {den_factored}")
    
    # Check if it's constant * (v - pole)
    if isinstance(den_factored, sp.Mul):
        args = den_factored.args
        const = 1
        for arg in args:
            if not sp.Symbol('v') in arg.free_symbols:
                const = const * arg
            else:
                # Should be (v - pole)
                pass
        print(f"  constant factor: {const}")
        
except Exception as e:
    print(f"Error: {e}")