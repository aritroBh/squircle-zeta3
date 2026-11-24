import sympy as sp
r5 = sp.sqrt(5)
rho = sp.sqrt(10-2*r5)
tm = 3*rho/4
print("as_powers_dict:", tm.as_powers_dict())
print("get rho:", tm.as_powers_dict().get(rho,0))
e = sp.expand(2*(rho/4) + 2*(rho*(1+r5)/8))
print("expanded:", e)
for t2 in sp.Add.make_args(e):
    print(" term:", t2, t2.as_powers_dict())
