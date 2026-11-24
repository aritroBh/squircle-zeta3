import sys, os, time
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
t0 = time.time()
# simplify each L(r) to Q(i) individually
Lc = {}
for r in range(1, 9):
    Lr = T.L(r)
    L2 = sp.radsimp(sp.simplify(Lr))
    assert not L2.has(sp.sqrt(2)), r
    Lc[r] = L2
    print(f'L({r}) simplified: ops={sp.count_ops(L2)} ({round(time.time()-t0,1)}s)', flush=True)
import pickle
pickle.dump({r: str(Lc[r]) for r in Lc}, open('L_simplified.pkl','wb'))
print('L cache saved')
