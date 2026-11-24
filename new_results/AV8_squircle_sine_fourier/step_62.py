import sympy as sp
j, k, m = sp.symbols('j k m', integer=True, positive=True)
# Symbolic evaluation of Kt(2k) = 2*sum (-1)^{j-1} K(2k-2j+1)
# K(1) = pi/8 - pi ln2/4 ; K(r) = -pi/(2(r^2-1)) for r>=3
# Write Kt(2k) = 2*(-1)^{k-1}*K(1) + 2*sum_{j=1}^{k-1} (-1)^{j-1} * (-pi/(2((2k-2j+1)^2-1)))
r = 2*k-2*j+1
term = 2*(-1)**(j-1) * (-sp.pi/(2*(r**2-1)))
S = sp.summation(term, (j, 1, k-1))
print("sum part:", sp.simplify(S))
