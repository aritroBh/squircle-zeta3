import sympy as sp, pickle, os, sys
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
x = sp.Symbol('x', positive=True)
z2, z3 = sp.Symbol('z2'), sp.Symbol('z3')
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)
def simp(e):
    e = e.replace(lambda f: isinstance(f, sp.Function) and f.func==sp.polylog and f.args[0]==1,
                  lambda f: -sp.log(1-f.args[1]))
    return sp.simplify(sp.expand(e))
print("== classical identity certificates (F' == 0) ==")
# Euler: Li2(x)+Li2(1-x) = z2 - log x log(1-x)
E1 = Li2(x)+Li2(1-x) - (z2 - sp.log(x)*sp.log(1-x))
print("Euler d/dx:", simp(sp.diff(E1,x)))
# Landen: Li2(x)+Li2(x/(x-1)) = -log^2(1-x)/2
L1 = Li2(x)+Li2(x/(x-1)) + sp.log(1-x)**2/2
print("Landen d/dx:", simp(sp.diff(L1,x)))
# inversion: Li2(-x)+Li2(-1/x) = -z2 - log^2(x)/2
I1 = Li2(-x)+Li2(-1/x) + z2 + sp.log(x)**2/2
print("inversion d/dx:", simp(sp.diff(I1,x)))
# Li3 three-term Landen: Li3(x)+Li3(1-x)+Li3(x/(x-1)) = z3 + z2 log(1-x) - (1/2) log^2(1-x) log x + (1/6) log^3(1-x)
LT = Li3(x)+Li3(1-x)+Li3(x/(x-1)) - (z3 + z2*sp.log(1-x) - sp.log(1-x)**2*sp.log(x)/2 + sp.log(1-x)**3/6)
dLT = sp.expand(sp.diff(LT,x))
# reduce with Landen & Euler for Li2
dLT = dLT.replace(lambda f: isinstance(f, sp.Function) and f.func==sp.polylog and f.args[0]==1,
                  lambda f: -sp.log(1-f.args[1]))
dLT = sp.expand(dLT)
# substitute Li2(x/(x-1)) = -Li2(x) - log^2(1-x)/2 and Li2(1-x) = z2 - Li2(x) - log x log(1-x)
dLT2 = dLT.subs({sp.polylog(2, x/(x-1)): -sp.polylog(2,x) - sp.log(1-x)**2/2,
                 sp.polylog(2, 1-x): z2 - sp.polylog(2,x) - sp.log(x)*sp.log(1-x)})
print("Li3-Landen d/dx after Li2 reductions:", sp.simplify(sp.expand(dLT2)))
print("(pins: Euler/Landen at x->0+, inversion at x=1 with Li2(-1)=-z2/2 from duplication, Li3-Landen at x->0+)")
