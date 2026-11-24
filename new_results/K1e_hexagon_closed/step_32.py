
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
A1,B,kap,q0,V02,L2,lA0,lA1,lB,rr = sp.symbols('A1 B kap q0 V02 L2 lA0 lA1 lB r')
l1,l2a,l3a,l3l,z2,z3 = sp.symbols('l1 l2a l3a l3l z2 z3')
atom_of = {'Li3(x)':l3a,'Li3(x/(x-1))':l3l,'Li2(x)*ln(1-x)':l2a*l1,'ln^3(1-x)':l1**3,
           'Li2(x)':l2a,'ln^2(1-x)':l1**2,'ln(1-x)':l1,'zeta3':z3,'zeta2*ln(1-x)':z2*l1,'zeta2':z2,'1':sp.Integer(1)}
E = sp.Integer(0)
for name,(num,den) in clean.items():
    E += num/L2**sp.Rational(5,2)*atom_of[name]
E = sp.expand(E.subs({kap: rr-1}))
poly = sp.Poly(E, lA0,lA1,lB)
print("log-monomials:", len(poly.terms()))
org = {}
for exps, coef in poly.terms():
    mu = lA0**exps[0]*lA1**exps[1]*lB**exps[2]
    c = sp.factor(sp.expand(coef)*L2**sp.Rational(5,2))/L2**sp.Rational(5,2)
    org[exps] = c
    print(f"lA0^{exps[0]} lA1^{exps[1]} lB^{exps[2]}: {len(sp.Add.make_args(sp.expand(coef)))} terms")
pickle.dump(org, open(f"{WORK}/Earc_by_logmon.pkl","wb"))
print(f"saved ({time.time()-t0:.0f}s)")
