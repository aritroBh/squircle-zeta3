import sympy as sp

tau, s = sp.symbols('tau s', positive=True)
u = tau*s
Nc = 9  # expand to this order in tau
lc = sp.series(sp.log(sp.cos(sp.pi/4+u)), u, 0, Nc+2).removeO().expand()
ls = sp.series(sp.log(sp.sin(sp.pi/4+u)), u, 0, Nc+2).removeO().expand()
A = sp.expand(lc/tau)   # q log cos, q=1/tau
B = sp.expand(ls/tau)
# log h = tau*log(exp(A)+exp(B));  use log-sum-exp: = tau*(A+B)/2 + tau*log(2 cosh((A-B)/2))
halfsum = sp.expand((A+B)/2)
halfdiff = sp.expand((A-B)/2)
logh = sp.expand( tau*halfsum + tau*sp.log(2) + tau*sp.log(sp.cosh(halfdiff)) )
logh_ser = sp.series(logh, tau, 0, Nc+1).removeO().expand()
# log m for s>0: log sin(pi/4+tau s) = ls
logm = sp.expand(ls)
G = sp.expand(logh_ser - logm)   # log(h/m)
G_ser = sp.series(G, tau, 0, Nc+1).removeO()
print("G = log(h/m) coefficients g_k(s):")
for k in range(1, Nc+1):
    gk = sp.expand(G_ser.coeff(tau, k))
    if gk != 0:
        print(f"g_{k}(s) =", gk)
