import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_budget.py").read()

old = """    x = 1/(sqrt(2)*y)
    M4num = sqrt(2)*gamma(mpf(11)/2)*x**(-mpf(11)/2)
    M4 = M4num / q['H_lo']      # ~ C y^4"""
new = """    x = (1/y-Delta)/sqrt(2)   # t<0 rate: e^{-u(1/y-Delta)/sqrt2}
    M4num = sqrt(2)*gamma(mpf(11)/2)*x**(-mpf(11)/2)
    M4 = M4num / q['H_lo']      # ~ C y^4 (1-yDelta)^{-11/2}"""
assert old in src
src = src.replace(old, new)

old = """    # delta0 for right tail rate
    yD2 = y*Delta
    delta0 = yD2*(2+yD2)/(1+yD2)**2
    lamR = bR*delta0"""
new = """    # right tail rate from Lambda' bound: lam'(t) <= num(Delta)/D(Delta) - b =: -rDelta
    yD2 = y*Delta
    numD = y*(2*bR/y - (1+sqrt(2)*y) - y*Delta)   # y[theta - (1+sqrt2 y) - y Delta], theta = pi/(sqrt2 y) = 2b/y... careful: y*theta = 2b
    numD = 2*bR - y*(1+sqrt(2)*y) - y*y*Delta
    DD = (1+yD2)*(1+yD2+sqrt(2)*y)
    rDelta = bR - numD/DD
    assert rDelta > 0
    lamR = rDelta"""
assert old in src
src = src.replace(old, new)
open("dcy2_budget.py","w").write(src)
from mpmath import mp
mp.dps = 50
import importlib, dcy2_budget
importlib.reload(dcy2_budget)
dcy2_budget.P.kap0 = 2.48
for N in [400, 2000, 10**4, 10**6, 10**9]:
    q = dcy2_budget.quant(N)
    print(f"N={N:9d}  ratio={float(q['E_total']/q['budget']):9.4f}  E_mom={float(q['E_mom']):.3e}  E_tailmom={float(q['E_tailmom']):.3e}  Tails={float(q['Tails']):.3e}")
