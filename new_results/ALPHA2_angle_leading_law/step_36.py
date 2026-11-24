import pickle
from fractions import Fraction as Fr
import sympy as sp
Fd = pickle.load(open('<project>/research/08_theory/ALPHA2_angle_leading_law/Fd_data.pkl','rb'))
sqc = Fr(1,2)**1  # squircle c=1/sqrt2 -> t = 1; add point t=1, F=1
t = sp.symbols('t')
def tval(c): return (1-c**2)/c**2
for d in range(2, 8):
    pts = [(sp.Rational(1), sp.Rational(1))]  # squircle
    for c, v in Fd[d].items():
        pts.append((sp.Rational(c.numerator, c.denominator), sp.Rational(v.numerator, v.denominator)))
    pts = [( (1-x**2)/x**2, y) for (x,y) in pts]
    # fit polynomial of minimal degree
    for deg in range(1, 8):
        coeffs = sp.symbols(f'a0:{deg+1}')
        poly = sum(coeffs[i]*t**i for i in range(deg+1))
        eqs = [sp.Eq(poly.subs(t, tp), yp) for tp, yp in pts]
        sol = sp.solve(eqs, coeffs, dict=True)
        if sol:
            sol = sol[0]
            p = sp.nsimplify(sp.expand(poly.subs(sol)))
            # verify on all points
            ok = all(sp.simplify(p.subs(t,tp)-yp)==0 for tp,yp in pts)
            if ok:
                print(f"d={d}: F_d(t) = {p}   (deg {deg}, {len(pts)} pts)")
                break
    else:
        print(f"d={d}: no poly fit up to deg 7; pts={pts}")
