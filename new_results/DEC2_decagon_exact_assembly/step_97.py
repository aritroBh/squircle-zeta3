import sympy as sp
y = sp.Symbol('y', positive=True)
z2, z3 = sp.Symbol('z2'), sp.Symbol('z3')
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)
def Ap(z):
    return (sp.polylog(2,z) + sp.log(1-z)**2/2)/z
# L1: A(-y)+A(-1/y) = -z2 log y - log^3 y/3 + Li3(-y) + 2 Li3(-1/y) - log y Li2(-y) + z3
R_L1 = -z2*sp.log(y) - sp.log(y)**3/3 + Li3(-y) + 2*Li3(-1/y) - sp.log(y)*Li2(-y) + z3
G = -Ap(-y) + Ap(-1/y)/y**2 - sp.diff(R_L1,y)
G = G.replace(lambda f: isinstance(f, sp.Function) and f.func==sp.polylog and f.args[0]==1,
              lambda f: -sp.log(1-f.args[1]))
# inversion: Li2(-y) = -z2 - log^2 y /2 - Li2(-1/y); log(1+1/y) = log(1+y) - log y
G2 = sp.expand(G.subs({sp.polylog(2,-y): -z2 - sp.log(y)**2/2 - sp.polylog(2,-1/y),
                       sp.log(1+1/y): sp.log(1+y) - sp.log(y)}))
print("L1 derivative residual:", sp.simplify(G2))
import mpmath as mp
mp.mp.dps = 60
z3v = mp.zeta(3)
# pin at y=1 needs A(-1) = -5 z3/8: check numerically (classical Euler sum)
s = mp.nsum(lambda n: (-1)**n*mp.harmonic(n)/n**2, [1, mp.inf])
print("A(-1) + 5 z3/8 (numeric):", mp.nstr(s + 5*z3v/8, 5))
