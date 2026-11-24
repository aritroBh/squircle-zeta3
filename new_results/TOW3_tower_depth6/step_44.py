work = '<project>/research/08_theory/TOW3_tower_depth6'
lib = r'''"""tow3lib.py -- TOW3 library: e-family closed forms at depths <= 6, Newton feed Lambda,
exact tower data interface. Everything rational in (d, m) over Q(i).
Status of ingredients:
  a_delta (delta<=4): PROVED in TOW2.  a_5, a_6: PROVED here (Bell from proved R_1..R_4), 720 cells each.
  e(2,1),e(3,1),e(4,1),e(4,2): PROVED in TOW2 (verified identities + data).
  e(5,2), e(6,3): zigzag tops (equation PROVED all s in TOW2; numerator verified vs data s<=6 there, and 104 keys here).
  e(5,1), e(6,2): TOW2 fits (CANDIDATES there); VERIFIED vs 414 fresh TOWN1 keys here, and
                  become PROVED once every identity they enter is verified (prove_V56, prove_adj56).
  e(6,1): fit here (414 keys, 0 bad); same upgrade path.
"""
import sympy as sp
import tow2lib as T
from collections import Counter

d, m = sp.symbols('d m')
i = sp.I

def Pi(l, dd=None):
    dd = dd if dd is not None else d
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*dd-2*jj-1)
    return out

# ---- a_delta (rho/e00) ----
_L = None
def _load_L():
    global _L
    if _L is None:
        import pickle, os
        _L = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'L_parsed.pkl'),'rb'))
    return _L
def _parts(rem, minpart, cur):
    if rem == 0:
        yield tuple(sorted(cur)); return
    for p in range(minpart, rem+1):
        yield from _parts(rem-p, p, cur+[p])
_a = {}
def arho(delta):
    if delta == 0: return sp.Integer(1)
    if delta in _a: return _a[delta]
    if delta <= 4:
        _a[delta] = T.arho_simple(delta)
        return _a[delta]
    L = _load_L()
    tot = sp.Integer(0)
    for part in _parts(delta, 1, []):
        C = Counter(part)
        coef = sp.Integer(1)
        for sz, mult in C.items():
            coef *= (L[sz]**mult) / sp.factorial(mult)
        tot += coef
    _a[delta] = sp.cancel(tot)
    return _a[delta]

# ---- zigzag tops ----
A = T.A
def e_top(delta):
    if delta % 2 == 0:
        s = delta//2
        return (-1)**s * A[s+1] / Pi(s)
    s = (delta-1)//2
    return (-1)**(s+1) * i * A[s+1] * (d-s-2) / Pi(s)

# ---- sub-top / level-1 closed forms ----
e21 = -1/Pi(1)
e31 = i*(d-3)/Pi(1)
e41 = (d**3 - 9*d**2 + sp.Rational(113,4)*d - i*d/2 + i*m - 29)/Pi(2)
e51 = (-4*i*d**4 + 52*i*d**3 + (-6-269*i)*d**2 + 12*d*m + (24+638*i)*d - 48*m - 552*i)/(12*Pi(2))
e61 = (-16*d**6 + 336*d**5 + (-3064+48*i)*d**4 - 96*i*d**3*m + (15424-576*i)*d**3 + 1152*i*d**2*m
       + (-44629+2364*i)*d**2 + (-48-4728*i)*d*m + (68999-3228*i)*d + 48*m**2 + 6456*i*m - 43398)/(96*Pi(3))
e62 = (-4*d**3 + 48*d**2 + (-199+2*i)*d - 4*i*m + 268)/(2*Pi(3))

_EF = {(2,1): e21, (3,1): e31, (4,1): e41, (5,1): e51, (6,1): e61, (6,2): e62}
def e(delta, s, dd=None, mm=None, conj=False):
    """e(delta,s)/e00 as rational function of (d,m); optional evaluation and i->-i."""
    if s > delta//2: return sp.Integer(0)
    if s == 0: val = arho(delta)
    elif s == delta//2 and not (delta,s) in _EF:
        val = e_top(delta)
    else:
        val = _EF[(delta,s)]
    if conj: val = val.subs({sp.I: -sp.I})
    if dd is not None: val = val.subs({d: dd, m: mm})
    return sp.expand(val)

# ---- tau constants (Q(i)), verified vs definitions in TOW2 ----
tau1 = T.tau1f; tau2 = T.tau2f

def rA0(dd=None, mm=None):
    dd = dd if dd is not None else d; mm = mm if mm is not None else m
    return i*(2*dd-1)/(mm+1)
def rB0(dd=None, mm=None):
    dd = dd if dd is not None else d; mm = mm if mm is not None else m
    return -i*(2*dd-1)/(dd-mm+1)

def Lam(delta, s, dd=None, mm=None, conj=False):
    """Newton feed for key (mm, dd-mm) off-diagonal parents; conj = i->-i throughout."""
    dd0 = dd if dd is not None else d
    mm0 = mm if mm is not None else m
    ii = -sp.I if conj else sp.I
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        t1 = tau1(delta-dp); t2 = tau2(delta-dp)
        if conj:
            t1 = (-1)**(delta-dp)*tau2(delta-dp); t2 = (-1)**(delta-dp)*tau1(delta-dp)
        cA = t1 * (mm0+1) / (2*dd0-dp-1)
        cB = t2 * (dd0-mm0+1) / (2*dd0-dp-1)
        rA = ii*(2*dd0-1)/(mm0+1); rB = -ii*(2*dd0-1)/(dd0-mm0+1)
        tot += cA * e(dp, s, dd0+1, mm0+1, conj) * rA
        tot += cB * e(dp, s, dd0+1, mm0, conj) * rB
    return sp.expand(tot)
'''
with open(f'{work}/tow3lib.py','w') as f: f.write(lib)
# quick self-test
import sys, os
sys.path.insert(0, work); os.chdir(work)
import tow3lib as T3
import sympy as sp
print('a5 =', sp.cancel(T3.arho(5)) == sp.sympify(__import__('pickle').load(open(f'{work}/arho56_simplified.pkl','rb'))[5], locals={'d':T3.d,'m':T3.m,'I':sp.I}))
print('e(5,2) =', T3.e(5,2))
print('e(6,3) =', T3.e(6,3))
print('Lam(5,2) sample at d=7,m=4:', sp.simplify(T3.Lam(5,2,7,4)))
print('check_tau:', __import__('tow2lib').check_tau(7))
