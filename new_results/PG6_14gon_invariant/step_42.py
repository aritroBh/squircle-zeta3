import sympy as sp, time
x = sp.symbols('x')
al_root = sp.RootOf(x**3+x**2-2*x-1, 0)  # which real root?
print("RootOf(,0) ~", sp.N(al_root, 20))
al_root2 = sp.RootOf(x**3+x**2-2*x-1, 2)
print("RootOf(,2) ~", sp.N(al_root2, 20))
# alpha = 2cos(2pi/7) ~ 1.24698 is the LARGEST; index 2 likely
al = sp.AlgebraicNumber(al_root2)
t0=time.time()
X1 = 1-al; X2 = 2-al-al**2; X5 = 2-al**2; X6 = 3-al-al**2
u = al**2+al-1
e = sp.expand((1-X1)**3 * (1-X2) + al**7)
print("expand test:", e, "time %.3f"%(time.time()-t0))
# test a division + expand
t0=time.time()
e2 = sp.expand( (X1/X6 + 1/(1-X1)) * (u**2) + X5/(1-X2) )
print("div test:", e2, "time %.3f"%(time.time()-t0))
# verify a known identity: 1-x1 = alpha  => 1-(1-alpha) = alpha OK; Landen x6=x1/(x1-1)
print("Landen check:", sp.simplify(X6 - X1/(X1-1)))
