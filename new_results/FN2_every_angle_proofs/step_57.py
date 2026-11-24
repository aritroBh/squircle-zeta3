with open("proof.md") as f:
    txt = f.read()
txt = txt.replace(
 "T_3=\\frac{d_2}{a}\\Bigl(2\\theta_{20}-\\theta_{11}+\\theta_{10}\\Bigr)+\\frac{d_3}{4a^2}\\theta_{10}",
 "T_3=\\frac{d_2}{a}\\Bigl(-2\\theta_{20}+\\theta_{11}-\\theta_{10}\\Bigr)+\\frac{d_3}{4a^2}\\theta_{10}")
with open("proof.md","w") as f:
    f.write(txt)
# verify the fix + independently re-verify T3's closed form with sympy
import sympy as sp
a, A, c, z = sp.symbols('a A c z')
th20 = -(a**2-A)/(2*a**3); th11 = -A/a**3; th10 = 1/a
d2 = 2*z*(1-z); d3 = -4*c*z*(1-z)
T3 = sp.simplify(d2/a*(-2*th20 + th11 - th10) + d3/(4*a**2)*th10)
print("T3 =", sp.simplify(T3 - (-4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3)))
print("fixed and re-verified")
