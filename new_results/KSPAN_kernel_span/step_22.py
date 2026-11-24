import sys, importlib, re, sympy as sp, fractions
sys.path.insert(0, "../WT2_weight2_all_k")

import fastred
importlib.reload(fastred)
from fastred import *

S2 = sp.sqrt(2)

def parse_rational_expr_fixed(expr_str):
    """Fixed parser that handles constant factors in denominator"""
    expr_str = expr_str.strip()
    
    try:
        expr = sp.sympify(expr_str)
        expr = sp.expand(expr)
        
        # Get numerator and denominator
        num, den = expr.as_numer_denom()
        
        # Factor denominator
        den_factored = sp.factor(den)
        
        # Extract constant factor and (v - pole) part
        const_factor = 1
        v_pole_part = None
        
        if isinstance(den_factored, sp.Mul):
            args = den_factored.args
            for arg in args:
                if not sp.Symbol('v') in arg.free_symbols:
                    const_factor = const_factor * arg
                else:
                    v_pole_part = arg
            if v_pole_part is None:
                # Could be just constant
                v_pole_part = 1
        else:
            v_pole_part = den_factored
        
        # Now solve v_pole_part = 0 for pole
        v = sp.Symbol('v')
        if v_pole_part == sp.Symbol('v'):
            pole_expr = 0
        else:
            pole_expr = sp.solve(sp.Eq(v_pole_part, 0), v)[0]
        
        pole_real = sp.re(pole_expr)
        pole_imag = sp.im(pole_expr)
        assert pole_imag == 0
        
        # Convert pole to Q(sqrt2)
        p0 = sp.simplify(pole_real.subs(S2, 0))
        p1 = sp.simplify((pole_real - p0)/S2)
        
        pole = QSI(fractions.Fraction(int(sp.Rational(p0).p), int(sp.Rational(p0).q)),
                  fractions.Fraction(int(sp.Rational(p1).p), int(sp.Rational(p1).q)))
        
        # Coefficient is num / const_factor
        coeff = sp.simplify(num / const_factor)
        
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

# Test
test_exprs = [
    "-1/(3*v + 3 + 3*sqrt(2))",
    "(8 - 7*sqrt(2))/(48*v)",
    "-7/(8*v + 8 + 8*sqrt(2))",
    "sqrt(2)/(16*v)"
]

for expr in test_exprs:
    print(f"\nParsing: {expr}")
    pole, coeff = parse_rational_expr_fixed(expr)
    print(f"  Pole: {pole}")
    print(f"  Coeff: {coeff}")