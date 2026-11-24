# Planner check of WT2r/WT2b: finite-difference identity, (I1), (I4), and the d = 2 values, exact.
import sympy as sp
from math import comb, factorial
I = sp.I
def dfact(n):  # double factorial, (-1)!! = 1
    return sp.factorial2(n) if n >= 0 else 1
ok = all(sum((-1)**s*comb(d, s)*(s**m if (s or m) else 1)*((d-s)**(d-m) if (d-s or d-m) else 1) for s in range(d+1)) == factorial(d)*(-1)**m
         for d in range(1, 40) for m in range(d+1))
print("finite-difference identity, d < 40, all m:", ok)
K = lambda d: (-1)**d*dfact(2*d-3)*I**(d+1)/(2*sp.sqrt(2)*sp.factorial(d))
print("(I1) K_{d+1}(d+1) = -i(2d-1)K_d, d<60:", all(sp.simplify(K(d+1)*(d+1) + I*(2*d-1)*K(d)) == 0 for d in range(2, 60)))
rho = lambda p, q: K(p+q)*(-1)**p*(1 + (1 if p > q else 0))*comb(p+q, p)
tau1 = (-1-I)/2
hhat = lambda p: (p+1)*tau1*rho(p+1, p)/(4*p-1)
print("(I4) Im hhat_p = Im rho_pp, p<60:", all(sp.simplify(sp.im(sp.expand(hhat(p))) - sp.im(sp.expand(rho(p, p)))) == 0 for p in range(1, 60)))
print("rho_20 =", sp.nsimplify(sp.expand(rho(2, 0))), " rho_11 =", sp.expand(rho(1, 1)), " hhat_1 =", sp.expand(hhat(1)))
