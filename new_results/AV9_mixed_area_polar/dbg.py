
import sympy as sp
import av9core as C
eta = C.eta
N = 3
# self-dual params through general_coeffs
alpha  = sp.series(-eta/(1+eta), eta, 0, N+1).removeO()
beta_  = eta
gamma  = sp.series(1/(1+eta), eta, 0, N+1).removeO()
P      = 2/(1+eta)
lnI, e = C.general_coeffs(N, alpha, beta_, gamma, P)
print("lnI[0] =", lnI[0])
print("e[0] =", e[0])
h, e2 = C.selfdual_coeffs(N)
for n in range(1, N+1):
    print(f"n={n}: general e[n] =", sp.expand(e[n]))
    print(f"      selfdual e[n]=", sp.expand(e2[n]))
    print("      ratio/diff:", sp.simplify(e[n]-e2[n]))
