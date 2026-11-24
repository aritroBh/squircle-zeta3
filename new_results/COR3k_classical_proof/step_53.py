import numpy as np
from cor3_common import *
import sympy as sp

print("=== Task 3: Li-argument ranges on (0,1) (need: in (-1,1] , monotone) ===")
grid = np.linspace(1e-9, 1 - 1e-9, 200001)
for nm in ARG_NAMES:
    f = sp.lambdify(x, ARGS[nm], 'numpy')
    v = f(grid)
    d = sp.diff(ARGS[nm], x)
    dv = sp.lambdify(x, d, 'numpy')(grid)
    mono = "increasing" if np.all(dv > 0) else ("decreasing" if np.all(dv < 0) else "NON-MONOTONE")
    print(f"  {nm:>14}: min = {v.min():.15f}, max = {v.max():.15f}, derivative sign: {mono}"
          f"  -> in (-1,1]: {v.min() > -1 and v.max() <= 1}")
