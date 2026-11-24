import sympy as sp, pickle, time
work = '<project>/research/08_theory/TOW5_tower_depth10'
import sys; sys.path.insert(0, work)
d, m = sp.symbols('d m')
i = sp.I
raw = pickle.load(open(work+'/R5_raw.pkl','rb'))
num5 = sp.sympify(raw['num5'], locals={'d':d,'m':m,'I':i})
den5 = sp.sympify(raw['den5'], locals={'d':d,'m':m,'I':i})
# Fully reduce R_5 = num5/den5 with cancel
R5 = sp.cancel(num5/den5)
nn, dd = sp.fraction(R5)
print("cancel'd R_5 denominator factored:", sp.factor(dd))
print("Pi_5 =", sp.factor(sp.prod(2*d-2*j-1 for j in range(1,6))))
print("dd / Pi_5 =", sp.simplify(dd/sp.prod(2*d-2*j-1 for j in range(1,6))))
print("numerator degree:", sp.degree(sp.expand(nn), d), "terms:", len(sp.Poly(sp.expand(nn), d, m).terms()))
# leading coeff of N_5
N5red = sp.expand(nn)
print("N_5 leading term (highest d):", sp.LT(sp.Poly(N5red, d, m), d) if False else "")
Pd = sp.Poly(N5red, d)
print("N_5 leading coeff in d:", sp.LC(Pd), " degree:", sp.degree(N5red, d))
