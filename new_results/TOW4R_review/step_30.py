
import os, sys, pickle, time
import sympy as sp
work = "<project>/research/08_theory/TOW4R_review"
sys.path.insert(0, work + "/tow4_copy")
import tow2lib as T2
d, m = T2.d, T2.m
raw = pickle.load(open(work + "/tow4_copy/L_simplified.pkl", 'rb'))
cache = {k: sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I}) for k, v in raw.items()}
print("cache keys:", sorted(cache.keys()))
t0 = time.time()
for r in range(1, 9):
    diff = sp.cancel(sp.together(T2.L(r) - cache[r]))
    print(f"L({r}): cache == tow2lib construction : {diff == 0}  ({round(time.time()-t0,1)}s)", flush=True)
# also compare L_simplified vs L_parsed (TOW3's cache)
raw2 = pickle.load(open(work + "/tow4_copy/L_parsed.pkl", 'rb'))
print("L_parsed keys:", sorted(raw2.keys()) if hasattr(raw2,'keys') else type(raw2))
