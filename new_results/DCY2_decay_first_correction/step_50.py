import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_budget.py").read()
old = """    # right tail rate from Lambda' bound: lam'(t) <= num(Delta)/D(Delta) - b =: -rDelta
    yD2 = y*Delta
    numD = y*(2*bR/y - (1+sqrt(2)*y) - y*Delta)   # y[theta - (1+sqrt2 y) - y Delta], theta = pi/(sqrt2 y) = 2b/y... careful: y*theta = 2b
    numD = 2*bR - y*(1+sqrt(2)*y) - y*y*Delta
    DD = (1+yD2)*(1+yD2+sqrt(2)*y)
    rDelta = bR - numD/DD
    assert rDelta > 0
    lamR = rDelta"""
new = """    # right tail rate: M(t) <= M(Delta) e^{-rDelta (t-Delta)}, rDelta = min(-Lambda'(Delta), b)
    yD2 = y*Delta
    uD_R = 1+yD2
    QR = uD_R**2 + sqrt(2)*y*uD_R + y**2
    LamPR = y*((N-1)/uD_R - N*(uD_R+y/sqrt(2))/QR) - bR
    assert LamPR < 0
    rDelta = min(-LamPR, bR)
    lamR = rDelta"""
assert old in src
src = src.replace(old, new)
old = """    def Pofu(u):
        return -2*bR*u**3 - 2*y*(1+pi/2)*u**2 + sqrt(2)*(N-2)*y**2*u + 2*(N-1)*y**3"""
new = """    def Pofu(u):
        return -2*bR*u**3 - (pi+2)*y*u**2 + (sqrt(2)*(N-2)*y**2 - 2*bR*y**3)*u + 2*(N-1)*y**3"""
assert old in src
src = src.replace(old, new)
open("dcy2_budget.py","w").write(src)
from mpmath import mp
mp.dps = 50
import importlib, dcy2_budget
importlib.reload(dcy2_budget)
dcy2_budget.P.kap0 = 2.48
for N in [400, 2000, 10**4, 15500, 10**5, 10**6, 10**9, 10**12]:
    q = dcy2_budget.quant(N)
    print(f"N={N:9d}  ratio={float(q['E_total']/q['budget']):9.4f}  E_mom={float(q['E_mom']):.3e}  E_tailmom={float(q['E_tailmom']):.3e}  Tails={float(q['Tails']):.3e}")
