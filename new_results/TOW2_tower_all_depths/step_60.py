import sys, os, time
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
t0 = time.time()
for delta in range(0,5):
    a = T.arho(delta)
    a2 = sp.radsimp(sp.simplify(a)) if delta >= 2 else a
    has2 = a2.has(sp.sqrt(2))
    print(f'arho({delta}): in Q(i): {not has2}  ({round(time.time()-t0,1)}s)', flush=True)
    if delta >= 2:
        # cache simplified
        pass
