# ======== STEP 1 indep. check #2: expand arclength integrand in u=e-1, integrate coefficients ========
import sympy as sp
u, lt, l1 = sp.symbols('u lt l1')
NN = 10
def T(poly, N):
    p = sp.expand(poly); return sum(p.coeff(u, n)*u**n for n in range(0, N+1))
# X = (1+u) exp(u lt),  Y = (1+u) exp(u l1)
ex = sp.exp(u*lt).series(u, 0, NN+1).removeO().expand()
ey = sp.exp(u*l1).series(u, 0, NN+1).removeO().expand()
X = T((1+u)*ex, NN); Y = T((1+u)*ey, NN)
W = T(X*X + Y*Y, NN)
# S = sqrt(W), recurrence S^2 = W
Sc = [sp.sqrt(W.coeff(u, 0))]
for n in range(1, NN+1):
    Sc.append(sp.expand((W.coeff(u, n) - sum(Sc[i]*Sc[n-i] for i in range(1, n)))/(2*Sc[0])))
S = sum(Sc[n]*u**n for n in range(NN+1))
open('step1_Sseries.txt','w').write(str([sp.expand(s) for s in Sc]))
print("s_n(t) = coefficient of u^n in arclength integrand sqrt(X^2+Y^2), as poly in lt=ln t, l1=ln(1-t):")
for n in range(0, 5): print("s_%d = %s" % (n, sp.expand(Sc[n])))
print("...") 
print("s_10 has", len(sp.Add.make_args(sp.expand(Sc[10]))), "terms")
