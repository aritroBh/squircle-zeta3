import pickle, sympy as sp
# Quick analysis: L_d leading coefficients and coefficient formulas for h_d
t = sp.symbols('t')
h = {2: t, 3: t+sp.Rational(1,3), 4: t**2+sp.Rational(2,3)*t, 5: t**2+t+sp.Rational(2,15),
     6: t**3+sp.Rational(4,3)*t**2+sp.Rational(17,45)*t, 7: t**3+sp.Rational(5,3)*t**2+sp.Rational(11,15)*t+sp.Rational(17,315)}
for d in sorted(h):
    m = d//2
    hh = sp.Poly(h[d], t)
    print(f"d={d}: h_d = {h[d]}   h(d,m-1)={hh.nth(m-1) if m-1>=0 else '-'} (formula (d-2)/3 = {sp.Rational(d-2,3)})")
