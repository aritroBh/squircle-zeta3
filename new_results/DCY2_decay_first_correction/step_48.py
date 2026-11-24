import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_budget.py").read()
old = """    def TailMom(p):
        # int_Delta^inf t^p e^{-c2 t^2 + |Re l1| t} dt with |Re l1| <= c2m Delta/2:
        lam = c2m*Delta/2
        # int_Delta^inf t^p e^{-lam t} dt = e^{-lam Delta} sum_{i=0}^p p!/(p-i)! Delta^{p-i} lam^{-i-1}
        return exp(-lam*Delta)*sum(factorial(p)/factorial(p-i)*Delta**(p-i)*lam**(-i-1) for i in range(p+1))
    E_tailmom = sum(cf*TailMom(p) for p,cf in pp.items())"""
new = """    from mpmath import gammainc
    assert ReL1 <= c2m*Delta/2, (N, ReL1, c2m*Delta/2)
    def TailMom(p):
        # int_Delta^inf t^p e^{-c2 t^2 + |Re l1| t} dt <= int_Delta^inf t^p e^{-c2 t^2/2} dt
        # = (1/2)(c2m/2)^{-(p+1)/2} Gamma((p+1)/2, c2m Delta^2/2)
        return mpf(1)/2*(c2m/2)**(-(p+1)/2)*gammainc((p+1)/2, c2m*Delta**2/2, regularized=False)
    E_tailmom = sum(cf*TailMom(p) for p,cf in pp.items())"""
assert old in src
src = src.replace(old, new)
open("dcy2_budget.py","w").write(src)
from mpmath import mp
mp.dps = 50
import importlib, dcy2_budget
importlib.reload(dcy2_budget)
dcy2_budget.P.kap0 = 2.48
for N in [400, 1000, 2000, 10**4, 15500, 10**5, 10**6, 10**9, 10**12]:
    q = dcy2_budget.quant(N)
    ratio = q['E_total']/q['budget']
    print(f"N={N:9d}  ratio={float(ratio):9.4f}   E_mom={float(q['E_mom']):.3e}  E_tailmom={float(q['E_tailmom']):.3e}  Tails={float(q['Tails']):.3e}")
