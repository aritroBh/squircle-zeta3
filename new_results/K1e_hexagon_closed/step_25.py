
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
simp = pickle.load(open(f"{WORK}/Earc_basis_A0sub.pkl","rb"))
A1,B,kap,q0,V02,Ls = sp.symbols('A1 B kap q0 V02 Ls')
Lp2 = Ls**2
def toLs(c):
    return sp.expand(c.subs({sp.Symbol('L2'): Ls**2}))
fact = {}
for name, c in simp.items():
    c2 = toLs(c)
    # common denominator
    num, den = sp.fraction(sp.together(c2))
    num = sp.factor(sp.expand(num))
    den = sp.factor(den)
    fact[name] = (num, den)
    pole = None
    dk = sp.Poly(den, kap).degree() if den.has(kap) and den.is_polynomial(kap) else 0
    print(f"{name:>16}: num #terms {len(sp.Add.make_args(sp.expand(num)))}, den = {den}   ({time.time()-t0:.0f}s)")
pickle.dump(fact, open(f"{WORK}/Earc_basis_fact.pkl","wb"))
print("saved")
