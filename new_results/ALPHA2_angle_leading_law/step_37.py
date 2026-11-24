import pickle, sympy as sp
from fractions import Fraction as Fr
Fd = pickle.load(open('<project>/research/08_theory/ALPHA2_angle_leading_law/Fd_data.pkl','rb'))
t = sp.symbols('t')
def tval(c): return (1-c**2)/c**2
fits = {}
for d in range(2, 8):
    pts = [(sp.Rational(1), sp.Rational(1))]
    for c, v in Fd[d].items():
        x = sp.Rational(c.numerator, c.denominator)
        pts.append(((1-x**2)/x**2, sp.Rational(v.numerator, v.denominator)))
    found = None
    for deg in range(0, 6):
        if len(pts) < deg + 2: continue  # need overdetermined
        coeffs = sp.symbols(f'a0:{deg+1}')
        poly = sum(coeffs[i]*t**i for i in range(deg+1))
        # solve on first deg+1 points, check the rest
        sol = sp.solve([sp.Eq(poly.subs(t,tp), yp) for tp,yp in pts[:deg+1]], coeffs, dict=True)
        if not sol: continue
        p = sp.expand(poly.subs(sol[0]))
        if all(sp.simplify(p.subs(t,tp)-yp)==0 for tp,yp in pts):
            found = p; break
    if found is not None:
        fits[d] = found
        print(f"d={d}: F_d(t) = {found}   verified on {len(pts)} pts (deg {deg})")
    else:
        print(f"d={d}: NO fit deg<=5, {len(pts)} pts")
pickle.dump(fits, open('<project>/research/08_theory/ALPHA2_angle_leading_law/Fd_fits.pkl','wb'))
