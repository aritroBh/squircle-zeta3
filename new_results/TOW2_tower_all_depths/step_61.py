import sys, os, time
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
t0 = time.time()
# fully simplify arho(1..4) to Q(i) and cache
cache = {}
for delta in range(1,5):
    a = T.arho(delta)
    a2 = sp.radsimp(sp.simplify(a))
    assert not a2.has(sp.sqrt(2)), delta
    cache[delta] = a2
    print(f'arho({delta}) simplified ({round(time.time()-t0,1)}s): {sp.collect(a2, [T.d, T.m]) if delta<=1 else "big"}', flush=True)
import pickle
pickle.dump({de: str(cache[de]) for de in cache}, open('arho_simplified.pkl','wb'))
print('saved')
# append cache loader to tow2lib
add = '''

# ---------------- simplified arho cache (Q(i) forms, derived from complete_bell) ----------------
import pickle as _pkl, os as _os
_ac = _pkl.load(open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'arho_simplified.pkl'), 'rb'))
_arho_simple = {}
def arho_simple(delta):
    if delta == 0: return sp.Integer(1)
    if delta not in _arho_simple:
        _arho_simple[delta] = sp.sympify(_ac[delta], locals={'d': d, 'm': m, 'I': sp.I})
    return _arho_simple[delta]
'''
txt = open('tow2lib.py').read()
open('tow2lib.py','w').write(txt + add)
print('patched tow2lib with arho_simple')
