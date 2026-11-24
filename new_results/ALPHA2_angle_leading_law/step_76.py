import sympy as sp
t, x = sp.symbols('t x')
f = {2: t, 3: 1+3*t, 4: 2*t+3*t**2, 5: 2+15*t+15*t**2,
     6: 17*t+60*t**2+45*t**3, 7: 17+231*t+525*t**2+315*t**3, 8: 62*t+378*t**2+630*t**3+315*t**4}
A = {2:1,3:2,4:5,5:16,6:61,7:272,8:1385}
print("P_d(x) (ascending coeffs), F_d = P_d/[A_d (1+ca)^floor(d/2)]:")
for d in range(2, 9):
    fd = f[d]; f1 = fd.subs(t, 1)
    P = sp.together(fd.subs(t, (1-x)/(1+x)))
    Pd = sp.together(P * (1+x)**(d//2))
    Pd = sp.expand(Pd * A[d] / f1)
    coeffs = sp.Poly(Pd, x).all_coeffs()[::-1]
    print(f"d={d}: {coeffs}")
