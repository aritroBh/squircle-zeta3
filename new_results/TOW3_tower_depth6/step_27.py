import sys, os, time, pickle
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
import sympy as sp
import tow2lib as T

t0 = time.time()
a5 = T.arho(5)
print('a5 assembled', round(time.time()-t0,1), 's', flush=True)
a6 = T.arho(6)
print('a6 assembled', round(time.time()-t0,1), 's', flush=True)
# simplify to single fractions
a5s = sp.cancel(sp.together(a5))
print('a5 simplified', round(time.time()-t0,1), 's', flush=True)
a6s = sp.cancel(sp.together(a6))
print('a6 simplified', round(time.time()-t0,1), 's', flush=True)
n5, d5 = sp.fraction(a5s); n6, d6 = sp.fraction(a6s)
print('a5 num deg', sp.Poly(n5, T.d).degree(), 'den:', sp.factor(d5))
print('a6 num deg', sp.Poly(n6, T.d).degree(), 'den:', sp.factor(d6))
pickle.dump({5: str(a5s), 6: str(a6s)}, open('arho56_simplified.pkl','wb'))
print('saved', round(time.time()-t0,1))
