import sys, os, time, pickle
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import sympy as sp
import tow2lib as T
work = '<project>/research/08_theory/TOW4_tower_depth8'
d, m = T.d, T.m
Ls_raw = pickle.load(open(work + '/L_simplified.pkl','rb'))
cache = {k: sp.sympify(v, locals={'d':d,'m':m,'I':sp.I}) for k,v in Ls_raw.items()}
allok = True
for r in range(1,9):
    t0=time.time()
    Lr = sp.cancel(T.L(r))
    diff = sp.cancel(sp.together(Lr - cache[r]))
    ok = (diff == 0)
    allok &= ok
    print(f"L({r}): cache == tow2lib.L proved construction : {ok}  ({round(time.time()-t0,1)}s)", flush=True)
print("ALL r=1..8:", allok)
