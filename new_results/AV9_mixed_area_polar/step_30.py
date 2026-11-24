
import sympy as sp
z2 = sp.zeta(2)
e = sp.pi**4/720 + sp.pi**2/2
print("subs pi**2:", e.subs(sp.pi**2, 6*z2))
print("expand coeff z2^2:", sp.expand((36*z2**2*sp.Rational(1,720))).coeff(z2**2))
# try dict loop
for pw, rep in {sp.pi**2: 6*z2, sp.pi**4: 36*z2**2}.items():
    e = e.subs(pw, rep)
print("after loop:", e)
