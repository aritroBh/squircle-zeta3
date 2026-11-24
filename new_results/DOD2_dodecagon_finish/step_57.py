import pickle, sympy as sp
import mpmath as mp, time
Q = pickle.load(open('basis_quotient.pkl','rb'))
maxden = 1
for pv, expr in Q['elim'].items():
    for j,c in expr.items():
        maxden = max(maxden, sp.Rational(c).q)
print("max elim denominator:", maxden)
maxnum = max(abs(sp.Rational(c).p) for pv,expr in Q['elim'].items() for c in expr.values())
print("max elim numerator:", maxnum)
