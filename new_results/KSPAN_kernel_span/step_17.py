import sys, importlib, re, sympy as sp, fractions
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *

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

# Test again
test_exprs = [
    "sqrt(2)/(16*v)",
    "-7/(8*v + 8 + 8*sqrt(2))",
    "(2*sqrt(2) + 7)/(16*v)",
    "-5/(16*v - 16 + 16*sqrt(2))"
]

for expr in test_exprs:
    print(f"\nParsing: {expr}")
    pole, coeff = parse_rational_expr(expr)
    print(f"  Pole: {pole}")
    print(f"  Coeff: {coeff}")