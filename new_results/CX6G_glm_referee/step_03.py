#!/usr/bin/env python3
"""Step 03: Verify display (3) - M(eps) integral forms"""
import mpmath as mp
mp.mp.dps = 80

def B(z):
    return 1 + mp.exp(-z) - mp.sqrt(1 + mp.exp(-2*z))

def C(z):
    return mp.exp(z) * B(z)

def M_eps_form1(eps):
    def integrand(x):
        return mp.exp(-eps*x) * (1+mp.exp(-x))**(-1-eps) * B((1-eps)*x)
    return -mp.quad(integrand, [0, mp.inf])

def M_eps_form2(eps):
    def integrand(x):
        return mp.exp(-x) * (1+mp.exp(-x))**(-1-eps) * C((1-eps)*x)
    return -mp.quad(integrand, [0, mp.inf])

print("Display (3) - M(eps) forms equivalence:")
for eps_val in [0.01, 0.05, 0.1, 0.2]:
    f1 = M_eps_form1(eps_val)
    f2 = M_eps_form2(eps_val)
    print(f"eps={eps_val}: diff={mp.nstr(abs(f1-f2), 10)}")
