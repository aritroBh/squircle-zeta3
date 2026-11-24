import pickle, sympy as sp
G2 = "<project>/research/08_theory/G2_M4alpha/"
n = sp.symbols('n')
terms = pickle.load(open(G2 + "mu4_formula.pkl", "rb"))
c4 = sp.Rational(0); c22 = sp.Rational(0)
for mon, cf in terms:
    e2, e3, e4 = mon[8], mon[9], mon[10]
    if e4 == 1: c4 += sp.sympify(cf)
    elif e2 == 2: c22 += sp.sympify(cf)
t1 = (n + 3)/((2*n + 3)*(4*n**2 - 1)); t2 = -n/((2*n + 3)*(4*n**2 - 1))
print("[z4]  mu identity:", sp.simplify(c4 - t1) == 0)
print("[z2^2] mu identity:", sp.simplify(c22 - t2) == 0)
# F1,F2 generating functions match these coefficients for all n by partial fractions (step_07); recheck apart decomposition
print("apart [z4]:", sp.apart(t1, n))
print("apart [z2^2]:", sp.apart(t2, n))
