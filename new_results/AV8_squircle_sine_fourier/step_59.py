import sympy as sp
m = sp.symbols('m', integer=True, positive=True)
t1 = 1/(64*m**2*(2*m+1)**4)
t2 = 1/(64*(m+1)**2*(2*m+1)**4)
# do partial fractions then sum piece by piece
pf1 = sp.apart(t1, m)
pf2 = sp.apart(t2, m)
print("pf1:", pf1)
print("pf2:", pf2)
