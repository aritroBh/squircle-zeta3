import sympy as sp
v=sp.symbols('v'); sq2=sp.sqrt(2)
Gm30_3 = (v**6 + 3*sq2*v**6 + 6*sq2*v**5 + 3*v**4 + 3*sq2*v**4 - 4*sq2*v**3 - 3*sq2*v**2 + 3*v**2 + 6*sq2*v - 3*sq2 + 1)/(24*(v**6+3*v**4+3*v**2+1))
Gm30_18 = (-29*v**6 + 3*sq2*v**6 + 66*sq2*v**5 - 87*v**4 + 3*sq2*v**4 + 116*sq2*v**3 - 87*v**2 - 3*sq2*v**2 + 66*sq2*v - 29 - 3*sq2)/(24*(v**6+3*v**4+3*v**2+1))
Gm30_20 = (-33*v**6 + 3*sq2*v**6 + 74*sq2*v**5 - 99*v**4 + 3*sq2*v**4 + 132*sq2*v**3 - 99*v**2 - 3*sq2*v**2 + 74*sq2*v - 33 - 3*sq2)/(24*(v**6+3*v**4+3*v**2+1))
s30 = -(sp.Rational(1,12))*(1 - 2*sq2*v/(v**2+1))
d18 = sp.simplify(Gm30_18 - Gm30_3)
d20 = sp.simplify(Gm30_20 - Gm30_3)
print("Gm30^18 - Gm30^3 =", sp.simplify(d18))
print("15*s30           =", sp.simplify(15*s30))
print("affine law k=18 exact:", sp.simplify(d18-15*s30)==0)
print("affine law k=20 exact:", sp.simplify(d20-17*s30)==0)
with open("<project>/research/08_theory/TOW1R_review/progress.log","a") as fh:
    fh.write("step_26: affine Gm law Gm30^(k)=Gm30^(3)+(k-3)s30 EXACT at k=18,20\n")
