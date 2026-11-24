
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
fact = pickle.load(open(f"{WORK}/Earc_basis_fact.pkl","rb"))
Ls = sp.Symbol('Ls')
for name in ['zeta3','Li3(x)','Li3(x/(x-1))','Li2(x)*ln(1-x)','zeta2*ln(1-x)','zeta2','Li2(x)']:
    num, den = fact[name]
    # check parity in Ls
    odd = sp.expand(num).subs({Ls:-Ls}) - sp.expand(num)
    print(f"=== {name}  (odd-in-Ls part zero: {sp.expand(odd)==0}) ===")
    print(sp.factor(num))
