import sympy as sp, time
t0 = time.time()
W, mu = sp.symbols('W mu')
a, A, B, D4 = sp.symbols('a A B D4')   # f', f'', f''', f'''' at 1
c, z = sp.symbols('c z')
k2, k3, k4, Xi = sp.symbols('k2 k3 k4 Xi')

# A_cal to order 4, Btilde to order 3
t = sp.Symbol('t')
logf = -sp.log(1 - a*t + A*t**2/2 - B*t**3/6 + D4*t**4/24)
Ac = sp.expand(sp.series(logf, t, 0, 5).removeO())
fp = a - A*t + B*t**2/2 - D4*t**3/6
bt = sp.expand(sp.series(sp.log(fp/a), t, 0, 4).removeO())

# Theta to total order 4
idx = [(i,j) for k in range(1,5) for i in range(1,k+1) for j in [k-i]]
vs = sp.symbols(' '.join(f'v{i}{j}' for i,j in idx))
Theta = sum(v*v**0 for v in [])  # build
Theta = 0
for (i,j), v in zip(idx, vs):
    Theta += v*W**i*mu**j
Acal = sum(sp.expand(Ac).coeff(t,p)*Theta**p for p in range(1,5))
Btil = sum(sp.expand(bt).coeff(t,p)*Theta**p for p in range(1,4))
eq = sp.expand(Acal - mu*Btil - W)
sols = {}
for (i,j), v in zip(idx, vs):
    coeff = sp.expand(eq.subs(sols)).coeff(W,i).coeff(mu,j)
    sols[v] = sp.solve(coeff, v)[0]
th = {ij: sols[v] for ij, v in zip(idx, vs)}

# Compare with P3's printed theta_{40}, theta_{31}, theta_{22}, theta_{13} (b=A, c=B, d=D4)
b_, c_, d_ = A, B, D4
th40_P3 = -(a**6 - 7*a**4*b_ - 6*a**3*c_ + 18*a**2*b_**2 - a**2*d_ + 10*a*b_*c_ - 15*b_**3)/(24*a**7)
th31_P3 = -(7*a**4*b_ + 6*a**3*c_ - 24*a**2*b_**2 + a**2*d_ - 13*a*b_*c_ + 23*b_**3)/(6*a**7)
th22_P3 = -3*b_*(2*a**2*b_ + a*c_ - 3*b_**2)/(2*a**7)
th13_P3 = -b_**3/a**7
print("theta40 match P3:", sp.simplify(th[(4,0)] - th40_P3) == 0)
print("theta31 match P3:", sp.simplify(th[(3,1)] - th31_P3) == 0)
print("theta22 match P3:", sp.simplify(th[(2,2)] - th22_P3) == 0)
print("theta13 match P3:", sp.simplify(th[(1,3)] - th13_P3) == 0)

# T4 from the master formula; tail d2=2z(1-z), d3=-4cz(1-z), d4=(1-c^2)(5c^2-1)*3/8
d2 = 2*z*(1-z); d3 = -4*c*z*(1-z); d4 = 3*(1-c**2)*(5*c**2-1)/8
# express d4 in z,c: 1-c^2 = 4z(1-z); c^2 = (1-2z)^2
T4 = ( (th[(4,0)]*sp.factorial(4)*(-1)**5*d2/a + th[(3,1)]*sp.factorial(3)*(-1)**4*d2/a
       + th[(2,2)]*sp.factorial(2)*(-1)**3*d2/a + th[(1,3)]*sp.factorial(1)*(-1)**2*d2/a)
     + (th[(3,0)]*sp.factorial(3)*(-1)**4*d2/a*(-1) + th[(2,1)]*sp.factorial(2)*(-1)**3*d2/a*(-1)
       + th[(1,2)]*sp.factorial(1)*(-1)**2*d2/a*(-1))
     + (th[(2,0)]*sp.factorial(2)*(-1)**3*d2/a*(1) + th[(1,1)]*sp.factorial(1)*(-1)**2*d2/a*(1))
     + (th[(1,0)]*sp.factorial(1)*(-1)**2*d2/a*(-1))
     # l=3: d3/(a^2 2^{i+1}), i+j+r=2, C(-2,r)=(-1)^r(r+1)
     + (th[(2,0)]*sp.factorial(2)*(-1)**3*d3/(a**2*2**3) + th[(1,1)]*sp.factorial(1)*(-1)**2*d3/(a**2*2**2))
     + (th[(1,0)]*sp.factorial(1)*(-1)**2*d3/(a**2*2**2)*(-2))
     # l=4: d4/(a^3 3^2)
     + th[(1,0)]*d4/(a**3*9) )
T4 = sp.expand(T4)
# Check c=0 (z=1/2) vs P3's T4 = (12 a^4 - 48 a^2 b - a^2 - 24 a c + 96 b^2)/(24 a^6)
T4_P3 = (12*a**4 - 48*a**2*A - a**2 - 24*a*B + 96*A**2)/(24*a**6)
print("T4 at c=0,z=1/2 vs P3:", sp.simplify(T4.subs({c:0, z:sp.Rational(1,2)}) - T4_P3) == 0)
print("T4 (general) =", sp.factor(T4))
print("elapsed", time.time()-t0)
