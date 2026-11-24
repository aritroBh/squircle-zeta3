import sys, importlib, re, sympy as sp, fractions
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *

# Helper to parse rational expressions from the reference output
def parse_rational_expr(expr_str):
    """Parse expression like 'sqrt(2)/(16*v)' or '(2*sqrt(2) + 7)/(16*v)'"""
    # Remove spaces
    expr_str = expr_str.strip()
    
    # Parse with sympy
    try:
        expr = sp.sympify(expr_str)
        
        # Convert to partial fraction form (simple poles only)
        # Poles are at: 0, a=1+sqrt(2), -a, b=sqrt(2)-1, -b
        # where a=AQ, -a=-AQ, b=BQ, -b=-BQ
        
        # We expect expr to be of form: coeff / (v - pole) or constant/v
        # Expand and check
        expr = sp.expand(expr)
        
        # Get numerator and denominator
        num, den = expr.as_numer_denom()
        
        # Check denominator format
        if den == sp.Symbol('v'):
            pole = ZERO
            coeff = num
        elif str(den).startswith('('):
            # Parse v +/- something
            # Pattern: (v +/- something) or (v + c) or (v - c)
            # where c could be a +/- b*sqrt(2)
            den_str = str(den)
            if den_str.startswith('(') and den_str.endswith(')'):
                den_str = den_str[1:-1]
            
            # Check for patterns like v + 1 + sqrt(2) or v - 1 + sqrt(2)
            # Actually should be v - pole
            # Try to solve v - pole = 0
            v = sp.Symbol('v')
            pole_expr = sp.solve(sp.Eq(den, 0), v)[0]
            
            # Convert pole to QSI format
            # pole could be: -1 - sqrt(2), 1 - sqrt(2), etc.
            pole_real = sp.re(pole_expr)
            pole_imag = sp.im(pole_expr)
            # Should be real
            assert pole_imag == 0
            
            # Convert to Q(sqrt2) number
            # pole_real is something like -1 - sqrt(2)
            # Extract rational + sqrt2*rational
            s2 = sp.sqrt(2)
            p0 = sp.simplify(pole_real.subs(s2, 0))
            p1 = sp.simplify((pole_real - p0)/s2)
            
            # Create QSI pole
            pole = QSI(fractions.Fraction(int(sp.Rational(p0).p), int(sp.Rational(p0).q)),
                      fractions.Fraction(int(sp.Rational(p1).p), int(sp.Rational(p1).q)))
            
            coeff = num
        else:
            # Try other patterns
            # Check if denominator is constant*v
            if str(den).endswith('*v'):
                # Pattern like 16*v
                coeff_const = sp.simplify(sp.S(str(den).replace('*v', '')))
                pole = ZERO
                coeff = sp.simplify(num / coeff_const)
            else:
                raise ValueError(f"Unhandled denominator: {den}")
        
        # Convert coeff to QSI
        coeff = sp.expand(coeff)
        c0 = sp.simplify(coeff.subs(s2, 0))
        c1 = sp.simplify((coeff - c0)/s2)
        
        coeff_qsi = QSI(fractions.Fraction(int(sp.Rational(c0).p), int(sp.Rational(c0).q)),
                       fractions.Fraction(int(sp.Rational(c1).p), int(sp.Rational(c1).q)))
        
        return pole, coeff_qsi
        
    except Exception as e:
        print(f"Error parsing '{expr_str}': {e}")
        raise

# Test with some examples
test_exprs = [
    "sqrt(2)/(16*v)",
    "-7/(8*v + 8 + 8*sqrt(2))",
    "(2*sqrt(2) + 7)/(16*v)"
]

for expr in test_exprs:
    print(f"\nParsing: {expr}")
    pole, coeff = parse_rational_expr(expr)
    print(f"  Pole: {pole}")
    print(f"  Coeff: {coeff}")