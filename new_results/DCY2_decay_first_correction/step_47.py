import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp
mp.dps = 50
import importlib, dcy2_budget
importlib.reload(dcy2_budget)
dcy2_budget.P.kap0 = 2.48
for N in [400, 2000, 10**4, 15500, 10**5, 10**6, 10**9, 10**12]:
    q = dcy2_budget.quant(N)
    ratio = q['E_total']/q['budget']
    print(f"N={N:9d}  ratio={float(ratio):9.4f}   E_mom={float(q['E_mom']):.3e}  E_tailmom={float(q['E_tailmom']):.3e}  Tails={float(q['Tails']):.3e}  V3={float(q['V3']):.3f}")
