import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from genintegrand import stirling2
from fractions import Fraction as Fr

# fit P_t(x) = S(x, x-t), polynomial of degree 2t, for t = 3, 4
def fit_poly(vals, deg):
    # vals: list of (x,y) Fractions; Vandermonde solve with Fractions
    n = deg+1
    A = [[Fr(x)**k for k in range(n)] for (x,y) in vals[:n]]
    b = [y for (x,y) in vals[:n]]
    # Gaussian elimination
    M = [row[:] + [b[i]] for i,row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c,n) if M[r][c] != 0)
        M[c],M[p] = M[p],M[c]
        M[c] = [v/M[c][c] for v in M[c]]
        for r in range(n):
            if r != c and M[r][c]:
                M[r] = [v - M[r][c]*w for v,w in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]

import sympy as sp
x = sp.Symbol('x')
P = {}
for t in range(0,5):
    vals = [(Fr(m), Fr(stirling2(m, m-t))) for m in range(0, 4*t+4)]
    coeffs = fit_poly(vals, 2*t)
    poly = sum(sp.Rational(c)*x**k for k,c in enumerate(coeffs))
    # verify against all values and the recurrence
    ok_vals = all(poly.subs(x, m) == stirling2(m, m-t) for m in range(0, 30))
    P[t] = sp.expand(poly)
    print(f'P_{t}(x) =', sp.factor(poly), ' vals ok:', ok_vals)
# verify recurrences symbolically
for t in range(1,5):
    rem = sp.expand(P[t].subs(x,x+1) - ((x+1-t)*P[t-1] + P[t]))
    print(f'recurrence t={t}: remainder', rem)
import pickle
pickle.dump({t: str(P[t]) for t in P}, open('stirling_tails.pkl','wb'))
