import sympy as sp
p = sp.Symbol('p')
i = sp.I
rem = -4*i*p**2 + 13*i*p + i*sp.binomial(4*p-6, 2)/2 - 21*i/2
rem2 = sp.expand(rem.replace(sp.binomial(4*p-6,2), (4*p-6)*(4*p-7)/2))
print('E2b-adjacent remainder after expanding binomial:', sp.simplify(rem2))
