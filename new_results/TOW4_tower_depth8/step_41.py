work = '<project>/research/08_theory/TOW4_tower_depth8'
lib = '''"""tow4lib.py -- TOW4 library: extends tow2lib/tow3lib to depths 7, 8.
Everything here is fit-free:
  - tau table extended to n = 8 (tau1(8) = -169/128 + 41 I/128, tau2(8) = conjugate), each
    verified EXACTLY against the tow2lib definitions tau1/tau2 (sqrt2 cancels; diff == 0) and
    numerically to 40 digits; conjugation symmetry conj(tau1_n) = (-1)^n tau2_n holds at n = 8.
  - arho(delta), delta <= 8: Bell from the L(1..8) cache, the cache PROVED equal to tow2lib.L(r)
    symbolically for r = 1..8 (step: L cache verification; R_1..R_4 proved in TOW2).
  - forced(delta, s): the FORCED ascending family F, F(delta,0) = a_delta,
    F(delta,s+1) = Lambda(delta,s) - F(delta,s); loads the cancel'd cache forced_F78.pkl built
    by forced78.py (closing conditions F(delta, floor(delta/2)+1) = 0 proved there for delta <= 8).
"""
import os, pickle
import sympy as sp
import tow2lib as T2

d, m = T2.d, T2.m
i = sp.I
_HERE = os.path.dirname(os.path.abspath(__file__))

# ---- tau to n = 8 ----
_TAU8 = {8: (sp.Rational(-169,128) + sp.Rational(41,128)*sp.I,
             sp.Rational(-169,128) - sp.Rational(41,128)*sp.I)}
def tau1f(n):
    if n in _TAU8: return _TAU8[n][0]
    return T2.tau1f(n)
def tau2f(n):
    if n in _TAU8: return _TAU8[n][1]
    return T2.tau2f(n)
def check_tau8():
    d1 = sp.simplify(sp.expand(T2.tau1(8) - _TAU8[8][0]))
    d2 = sp.simplify(sp.expand(T2.tau2(8) - _TAU8[8][1]))
    return d1 == 0 and d2 == 0

def Pi(l, dd=None):
    return T2.Pi(l, dd if dd is not None else d)

# ---- a_delta (proved), delta <= 8 ----
from collections import Counter
_Lraw = pickle.load(open(os.path.join(_HERE, 'L_simplified.pkl'), 'rb'))
_L = {k: sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I}) for k, v in _Lraw.items()}
def _parts(rem, minpart, cur):
    if rem == 0:
        yield tuple(sorted(cur)); return
    for p in range(minpart, rem + 1):
        yield from _parts(rem - p, p, cur + [p])
_a = {}
def arho(delta):
    if delta == 0: return sp.Integer(1)
    if delta in _a: return _a[delta]
    tot = sp.Integer(0)
    for part in _parts(delta, 1, []):
        C = Counter(part)
        coef = sp.Integer(1)
        for sz, mult in C.items():
            coef *= (_L[sz]**mult) / sp.factorial(mult)
        tot += coef
    _a[delta] = sp.cancel(tot)
    return _a[delta]

# ---- forced family (from forced78.py cache) ----
_Fcache = None
def _load_F():
    global _Fcache
    if _Fcache is None:
        raw = pickle.load(open(os.path.join(_HERE, 'forced_F78.pkl'), 'rb'))
        _Fcache = {k: sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I}) for k, v in raw.items()}
    return _Fcache
def forced(delta, s, dd=None, mm=None, conj=False):
    if s > delta // 2: return sp.Integer(0)
    F = _load_F()
    val = F[(delta, s)]
    if conj: val = val.subs({sp.I: -sp.I})
    if dd is not None: val = val.subs({d: dd, m: mm})
    return sp.expand(val)
'''
open(work + '/tow4lib.py','w').write(lib)
# test import and tau8 check
import sys
sys.path.insert(0, work)
import tow4lib as T4
print("check_tau8:", T4.check_tau8())
print("arho(7) matches a7_closed:", end=' ')
import pickle
a7 = sp.sympify(pickle.load(open(work+'/a7_closed.pkl','rb')), locals={'d':T4.d,'m':T4.m,'I':sp.I})
import sympy as sp
print(sp.cancel(T4.arho(7)-a7)==0)
a8 = sp.sympify(pickle.load(open(work+'/a8_closed.pkl','rb')), locals={'d':T4.d,'m':T4.m,'I':sp.I})
print("arho(8) matches a8_closed:", sp.cancel(T4.arho(8)-a8)==0)
