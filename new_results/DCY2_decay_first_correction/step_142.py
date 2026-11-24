import sympy as sp
u = sp.Symbol('u')
v = sp.sqrt((1-sp.exp(-2*u))/(2*u))
ser = sp.series(v, u, 0, 12).removeO().expand()
print(ser)
c = [ser.coeff(u, k) for k in range(12)]
for k, ck in enumerate(c):
    print(f"c_{k} = {ck} = {sp.N(ck, 10) if ck != 0 else 0}")
