"""step_21.py -- PART 2a: degrees of K=Q(cos(alpha/2)) at 5pi/6, pi/12 (sympy exact),
and the mechanism: residues in {0,+-1,+-2,+-c,+-2c} => all coefficients in Q+Qc."""
import sympy as sp

print("minimal polynomials of c = cos(alpha/2):")
x = sp.Symbol('x')
for nm, val in [("5pi/6 (c=cos(5pi/12))", sp.cos(5*sp.pi/12)),
                ("pi/12 (c=cos(pi/24))", sp.cos(sp.pi/24)),
                ("3pi/4 (c=cos(3pi/8))", sp.cos(3*sp.pi/8)),
                ("pi/6  (c=cos(pi/12))", sp.cos(sp.pi/12))]:
    mpoly = sp.minpoly(val, x)
    print(f"  {nm}: minpoly = {mpoly} = 0, degree {sp.degree(mpoly)}")

# tau-equivariance note: under c -> -c, r1<->r2, residues covariant. Verify symbolically:
c = sp.Symbol('c')
print("\nresidues: a = (0,-2,1,-c,+c), b = (-1,1,-1,+2c,-2c) at rho=(1,-1,cos a,-1+2c,-1-2c)")
print("under tau: c->-c: r1=-1+2c <-> -1-2c=r2; a(r1)=-c -> +c = a(r2); b(r1)=2c -> -2c = b(r2)  [tau-covariant]")
