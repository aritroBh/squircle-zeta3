tow5lib = '''"""tow5lib.py -- TOW5 library: extends tow2lib/tow3lib/tow4lib to depths 9, 10.
Everything fit-free:
  - R_5 = S^(5)(t*)/S(t*) PROVED (prove_S5.py): R5 = num5/den5 (R5_raw.pkl); sl_5 = (log S)^(5)(t*)
    via exact cumulant algebra from the PROVED R_1..R_5.
  - tau table extended to n = 10 (tau910.pkl), each verified EXACTLY vs tow2lib definitions
    (sqrt2-aware reduction, difference 0) + 40-digit numeric; conj(tau1_n) = (-1)^n tau2_n holds.
  - L cache extended to r = 9, 10 (L910_candidate.pkl): simplified Q(i)(d,m) forms, each EQUAL to
    the tow2lib-style L(r) construction (with sl_5) symbolically (sqrt2-aware, difference 0).
  - arho(delta), delta <= 10: Bell from the L cache.
  - forced(delta, s), delta <= 10: FORCED ascending family from forced910.py cache (forced_F910.pkl).
"""
import os, pickle
import sympy as sp
import tow2lib as T2
import tow4lib as T4

d, m = T2.d, T2.m
i = sp.I
_HERE = os.path.dirname(os.path.abspath(__file__))

# ---- R_5 (PROVED, prove_S5.py) ----
_raw5 = pickle.load(open(os.path.join(_HERE, 'R5_raw.pkl'), 'rb'))
_num5 = sp.sympify(_raw5['num5'], locals={'d': d, 'm': m, 'I': sp.I})
_den5 = sp.sympify(_raw5['den5'], locals={'d': d, 'm': m, 'I': sp.I})
R5 = sp.cancel(_num5 / _den5)

# sl_5 = (log S)^(5)(t*) by exact cumulant algebra from PROVED R_1..R_5
R = T2.R
sl5 = sp.cancel(24*R[1]**5 - 60*R[1]**3*R[2] + 20*R[1]**2*R[3] + 30*R[1]*R[2]**2
                - 5*R[1]*R[4] - 10*R[2]*R[3] + R5)
T2.sl[5] = sl5

def Pi(l, dd=None):
    return T2.Pi(l, dd if dd is not None else d)

# ---- tau to n = 10 ----
_tau910 = pickle.load(open(os.path.join(_HERE, 'tau910.pkl'), 'rb'))
_TAU = {n: tuple(sp.sympify(v, locals={'I': sp.I}) for v in pair) for n, pair in _tau910.items()}
def tau1f(n):
    if n in _TAU: return _TAU[n][0]
    return T4.tau1f(n)
def tau2f(n):
    if n in _TAU: return _TAU[n][1]
    return T4.tau2f(n)

# ---- L cache (r = 1..10) ----
_Lraw = pickle.load(open(os.path.join(_HERE, 'L_simplified.pkl'), 'rb'))
_L = {k: sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I}) for k, v in _Lraw.items()}
_L910 = pickle.load(open(os.path.join(_HERE, 'L910_candidate.pkl'), 'rb'))
for k, v in _L910.items():
    _L[k] = sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I})

from collections import Counter
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

# ---- forced family (delta <= 10, from forced910.py cache) ----
_Fcache = None
def _load_F():
    global _Fcache
    if _Fcache is None:
        raw = pickle.load(open(os.path.join(_HERE, 'forced_F910.pkl'), 'rb'))
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
work = '<project>/research/08_theory/TOW5_tower_depth10'
with open(work+'/tow5lib.py','w') as f:
    f.write(tow5lib)
print("tow5lib.py written")

# Now compute a_9, a_10
import subprocess, time
PY = '<project>/.venv/bin/python'
code = '''
import sys, time, pickle
sys.path.insert(0, %r)
import sympy as sp
import tow5lib as T5
d, m = T5.d, T5.m
t0=time.time()
for delta in [9, 10]:
    t1=time.time()
    a = T5.arho(delta)
    nn, dd = sp.fraction(sp.together(a))
    print(f"a_{delta}: built ({round(time.time()-t1,1)}s); den factored: {sp.factor(dd)}", flush=True)
    print(f"   numerator degree in d: {sp.degree(sp.expand(nn), d)}, terms: {len(sp.Poly(sp.expand(nn), d, m).terms())}", flush=True)
    pickle.dump(str(sp.expand(a)), open(%r + f"/a{delta}_closed.pkl","wb"))
print("total", round(time.time()-t0,1))
''' % (work, work)
r = subprocess.run(['nice','-n','10',PY,'-c',code], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-1500:] if r.stderr else "")
