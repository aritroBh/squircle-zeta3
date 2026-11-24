import sys, os, time
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
t0 = time.time()
cache = {}
import pickle
ac = pickle.load(open('arho_simplified.pkl','rb'))
for delta in (5,6):
    a = T.arho(delta)
    print(f'arho({delta}) built ({round(time.time()-t0,1)}s), simplifying...', flush=True)
    a2 = sp.radsimp(sp.simplify(a))
    assert not a2.has(sp.sqrt(2)), delta
    ac[delta] = str(a2)
    print(f'arho({delta}) simplified ({round(time.time()-t0,1)}s), ops={sp.count_ops(a2)}', flush=True)
pickle.dump(ac, open('arho_simplified.pkl','wb'))
print('cache saved')
