src = open("cor6_nf.py").read()
# replace the hand-rolled factor_minpoly_modp with sympy-based
idx = src.index("def factor_minpoly_modp")
src = src[:idx] + '''def primes_above(m, p):
    """number of distinct primes above p in Q(theta) (p not dividing index), and factor degrees."""
    import sympy as sp
    x = sp.symbols('x')
    poly = sp.Poly(*[int(c) for c in reversed(m)], x, modulus=p)
    fac = sp.polys.factories.factor_list(poly, modulus=p)
    # fac: (lc, [(f, mult), ...])
    degs = [(sp.degree(f[0]), f[1]) for f in fac[1]]
    return sum(mult for deg, mult in degs), degs
'''
open("cor6_nf.py", "w").write(src)

from fractions import Fraction as Fr
import importlib, cor6_nf
importlib.reload(cor6_nf)
from cor6_nf import NF, primes_above

# define targets: (label, minpoly ascending, c0 as polynomial in theta)
targets = {
    "c=1/2 (hexagon)":   ([-2, 0, 1], (Fr(1,2), Fr(0))),       # theta = sqrt2? no - use Q: handle separately
}
# Let's define each with explicit minpoly and c0 poly
def show_target(label, minpoly, c0):
    F = NF(minpoly)
    c = c0
    one = F.one()
    two = F.from_int(2)
    N2 = F.norm(two); Nc = F.norm(c)
    N1mc = F.norm(F.sub(one, c)); N1pc = F.norm(F.add(one, c))
    print(f"{label}: deg {F.d}, N(c)={Nc}, N(1-c)={N1mc}, N(1+c)={N1pc}")

print("quadratic/polygon corners:")
show_target("c=sqrt2/2 (square,n=4)", [-2, 0, 1], (Fr(0), Fr(1, 2)))
show_target("c=sqrt3/2 (triangle,n=3)", [-3, 0, 1], (Fr(0), Fr(1, 2)))
show_target("c=sin(pi/10) (n=10)", [-5, 0, 1], (Fr(-1, 4), Fr(1, 4)))
show_target("c=sin(pi/5) (n=5)", [5, 0, -20, 0, 16], (Fr(0), Fr(1)))
show_target("c=sin(pi/8)=cos(3pi/8) (n=8)", [1, 0, -8, 0, 8], (Fr(0), Fr(1)))
show_target("c=sin(pi/12) (n=12)", [1, 0, -16, 0, 16], (Fr(0), Fr(1)))
