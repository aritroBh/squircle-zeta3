work = '<project>/research/08_theory/TOW4_tower_depth8'
script = r'''"""forced78.py -- FORCED ascending recursion to depth 8, symbolic in (d,m), fit-free.
F(delta,0) := a_delta (proved: Bell from L(1..8), L cache == tow2lib.L symbolically r<=8).
F(delta,s+1) := Lambda(delta,s) - F(delta,s).
Closing condition: F(delta, floor(delta/2)+1) == 0  (the tower law content at depth delta, non-adjacent).
Saves F(delta,s) (cancel'd) to forced_F78.pkl incrementally; prints closing verdicts.
"""
import sys, os, time, pickle
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/TOW4_tower_depth8'
sys.path.insert(0, work)
import sympy as sp
import tow2lib as T2
d, m = T2.d, T2.m
I = sp.I
tau1 = T2.tau1f; tau2 = T2.tau2f
from collections import Counter

# ---- a_delta (proved) ----
Ls_raw = pickle.load(open(work + '/L_simplified.pkl','rb'))
L = {k: sp.sympify(v, locals={'d':d,'m':m,'I':sp.I}) for k,v in Ls_raw.items()}
def parts(rem, minpart, cur):
    if rem == 0:
        yield tuple(sorted(cur)); return
    for p in range(minpart, rem+1):
        yield from parts(rem-p, p, cur+[p])
_a = {}
def arho(delta):
    if delta == 0: return sp.Integer(1)
    if delta in _a: return _a[delta]
    tot = sp.Integer(0)
    for part in parts(delta,1,[]):
        C = Counter(part)
        coef = sp.Integer(1)
        for sz,mult in C.items():
            coef *= (L[sz]**mult)/sp.factorial(mult)
        tot += coef
    _a[delta] = sp.cancel(tot)
    return _a[delta]

Fz = {}
def F(dp, s):
    return Fz[(dp, s)]
def Lam(delta, s):
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        fA = F(dp, s).subs({d: d+1, m: m+1}, simultaneous=True)
        fB = F(dp, s).subs({d: d+1, m: m}, simultaneous=True)
        tot += I*(2*d-1)/(2*d-dp-1) * (tau1(delta-dp)*fA - tau2(delta-dp)*fB)
    return tot

def save():
    pickle.dump({k: str(v) for k,v in Fz.items()}, open(work + '/forced_F78.pkl','wb'))

t0 = time.time()
closing = {}
for delta in range(0, 9):
    Fz[(delta, 0)] = sp.cancel(arho(delta))
    for s in range(0, delta//2+1):
        ts = time.time()
        nxt = sp.together(Lam(delta, s) - F(delta, s))
        num = sp.expand(sp.numer(nxt))
        is0 = (num == 0)
        if s+1 <= delta//2:
            Fz[(delta, s+1)] = sp.cancel(nxt)
            print(f'  F({delta},{s+1}) built  ({round(time.time()-ts,1)}s, total {round(time.time()-t0,1)}s)', flush=True)
        else:
            closing[delta] = is0
            print(f'  CLOSING delta={delta}: F({delta},{s+1}) numerator == 0 : {is0}  ({round(time.time()-ts,1)}s, total {round(time.time()-t0,1)}s)', flush=True)
        save()
    print(f'delta={delta} done ({round(time.time()-t0,1)}s)', flush=True)
print('ALL closing conditions delta<=8:', closing, flush=True)
print('VERDICT non-adjacent depth<=8:', all(closing.get(dd, True) for dd in range(0,9)), flush=True)
'''
open(work + '/forced78.py','w').write(script)
print("wrote forced78.py", len(script))
