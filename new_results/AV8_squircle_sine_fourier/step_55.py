import sympy as sp
n = 497323236409786642155382248146820840100456150797347717440463976893159497012533375533056
fac = sp.factorint(n)
print(fac)
# so log(n) = sum p_i^e_i -> e_i ln p_i
ln_n = sum(e*sp.log(q) for q,e in fac.items())
polyg = sp.pi**2/2 - 4  # polygamma(1,3/2)
N2 = (-ln_n - 7*sp.pi**2 + 96*sp.log(2)**2 + 54*polyg + 216)*sp.Rational(1,2)/(48)
# beta(3/2,1/2)=pi/2, so beta/(48 pi) = 1/96
N2 = sp.simplify((-ln_n - 7*sp.pi**2 + 96*sp.log(2)**2 + 54*polyg + 216)/96)
print("N''(2) =", N2)
print(sp.N(N2, 30))
