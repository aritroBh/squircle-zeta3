import sympy as sp
x = sp.Symbol('x', positive=True)
z2, z3 = sp.Symbol('z2'), sp.Symbol('z3')
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)
# A'(x) = (Li2(x) + log^2(1-x)/2)/x  (proved atom (1,0,1))
def Ap(z):
    return (sp.polylog(2,z) + sp.log(1-z)**2/2)/z
# DIFF: A(x) - A(1-x) = 2Li3(x)-2Li3(1-x)+log(1-x)Li2(1-x)-log x Li2(x)-1/2 log^2 x log(1-x)+1/2 log^2(1-x) log x
RHS = 2*Li3(x)-2*Li3(1-x)+sp.log(1-x)*Li2(1-x)-sp.log(x)*Li2(x)-sp.log(x)**2*sp.log(1-x)/2+sp.log(1-x)**2*sp.log(x)/2
F = Ap(x) + Ap(1-x) - sp.diff(RHS, x)   # derivative of [A(x)-A(1-x)-RHS]
F = F.replace(lambda f: isinstance(f, sp.Function) and f.func==sp.polylog and f.args[0]==1,
              lambda f: -sp.log(1-f.args[1]))
# reduce Li2(1-x) via Euler: Li2(1-x) = z2 - Li2(x) - log x log(1-x)
F2 = sp.expand(F.subs(sp.polylog(2,1-x), z2 - sp.polylog(2,x) - sp.log(x)*sp.log(1-x)))
print("DIFF derivative residual:", sp.simplify(F2))
# pin at x=1/2: both sides 0 trivially (A(x)-A(1-x) and RHS antisymmetric) -> exact
print("pin at x=1/2: LHS-RHS antisymmetry -> 0 exactly")
# L2: A(x)+A(x/(x-1)) = z2 log(1-x) + 2 log(1-x) Li2(1-x) + (1/2) log^2(1-x) log x - 3 Li3(1-x) + 3 z3
R_L2 = z2*sp.log(1-x) + 2*sp.log(1-x)*Li2(1-x) + sp.log(1-x)**2*sp.log(x)/2 - 3*Li3(1-x) + 3*z3
u = x/(x-1)
G = Ap(x) + Ap(u)*sp.diff(u,x) - sp.diff(R_L2,x)
G = G.replace(lambda f: isinstance(f, sp.Function) and f.func==sp.polylog and f.args[0]==1,
              lambda f: -sp.log(1-f.args[1]))
G2 = sp.expand(G.subs({sp.polylog(2, x/(x-1)): -sp.polylog(2,x) - sp.log(1-x)**2/2,
                       sp.polylog(2, 1-x): z2 - sp.polylog(2,x) - sp.log(x)*sp.log(1-x),
                       sp.log(1 - x/(x-1)): -sp.log(1-x)}))
print("L2 derivative residual:", sp.simplify(G2))
# pin x->0+: A(0)+A(0)=0, R_L2 -> z2*0 + 0 + 0 - 3 z3 + 3 z3 = 0  exact
