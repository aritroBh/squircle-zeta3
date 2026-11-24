import sympy as sp
xi, eta, u, s = sp.symbols('xi eta u s')
# beta_w := u^w - (xi^w + eta^w) reduced modulo eta = u - xi, expressed in (u, s=xi*eta)
# Lucas: xi^w+eta^w = p_w(u,s), p_w = u p_{w-1} - s p_{w-2}
p = {0: sp.Integer(2), 1: u}
for w in range(2, 9):
    p[w] = sp.expand(u*p[w-1] - s*p[w-2])
for w in range(2, 9):
    beta = sp.expand(u**w - p[w])
    C_w = (-1)**(w+1)*sp.factorial(w-1)
    P_w = sp.expand(C_w/sp.factorial(w)*beta)
    print(f'beta_{w} = {sp.factor(beta)}')
    print(f'P_{w} = {P_w}')
