common = r'''"""OCT9c common: load OCT9b stored data, evaluate brackets and reduction combos at high precision."""
import os, sys, pickle, itertools
BASE = "<project>/research/08_theory"
O9B = os.path.join(BASE, "OCT9b_octagon_debug")
if O9B not in sys.path: sys.path.insert(0, O9B)
import sympy as sp, mpmath as mp
from edge_engine import EdgeEngine
from syminfra import S_Li3, S_mLi2logp, S_logprod, SymbolSpace, symbol_of_expr

def split_weight(expr, sym2key):
    """(from OCT9b step_84): z2,z3 count weight 0 -> land in Wlo."""
    W3, Wlo = sp.Integer(0), sp.Integer(0)
    for term in sp.Add.make_args(sp.expand(expr)):
        coeff, rest = term.as_coeff_Mul()
        w = 0
        for f in sp.Mul.make_args(rest):
            b, e = f.as_base_exp()
            key = sym2key[b]
            if key[0] in ('z2','z3'): w += 0
            else: w += {'log':1,'Li2':2,'Li3':3,'Li21':3}[key[0]]*int(e)
        if w == 3: W3 += term
        else: Wlo += term
    return sp.expand(W3), sp.expand(Wlo)

K12 = {0: (sp.Rational(3,2), sp.Rational(-7,12)), 1: (sp.Rational(1,3), sp.Rational(1,9)),
       2: (sp.Rational(5,3), sp.Rational(-5,9)), 3: (sp.Rational(1,4), sp.Rational(1,24))}
LEN2 = {0: sp.Rational(41,400), 1: sp.Rational(25,1296), 2: sp.Rational(325,5184), 3: sp.Rational(29,1600)}  # |y'|^2

def get_engine(e):
    k1,k2 = K12[e]
    return EdgeEngine(k1,k2, tag=f"o9c_e{e}")

# ---- load stored data ----
def load_all():
    D = {}
    D['brackets'] = pickle.load(open(O9B+'/brackets.pkl','rb'))
    D['s123'] = pickle.load(open(O9B+'/sol_B123.pkl','rb'))
    D['s0'] = pickle.load(open(O9B+'/sol_B0.pkl','rb'))
    D['lists'] = pickle.load(open(O9B+'/Sunit_lists.pkl','rb'))   # Gargs(99), Dargs(27) over {2,3,5}
    D['sols2357'] = pickle.load(open(O9B+'/sols2357.pkl','rb'))
    D['lam'] = pickle.load(open(O9B+'/lam_cache.pkl','rb'))       # 340-digit strings
    D['Wp0'] = pickle.load(open(O9B+'/W_parts_e0.pkl','rb'))      # W3_e0, Wlo_e0
    D['P3'] = pickle.load(open(O9B+'/P3_assembled.pkl','rb'))
    return D

# ---- Li21 label map (a,c)->cache key, from OCT9b step_62 ----
LAMKEY = {}
for e, lst in {0:[((-3,2),(-7,18)),((7,12),(-18,7)),((3,5),(-7,3)),((-7,5),(-3,7))],
               1:[((-1,3),(1,3)),((-1,9),(3,1)),((1,4),(2,5)),((1,10),(5,2))],
               2:[((-5,3),(-1,3)),((5,9),(-3,1)),((5,8),(-2,1)),((-5,4),(-1,2))],
               3:[((-1,4),(1,6)),((-1,24),(6,1)),((1,5),(5,29)),((1,29),(29,5))]}.items():
    for j,(a,c) in enumerate(lst):
        LAMKEY[(e, sp.Rational(*a), sp.Rational(*c))] = f"e{e}L{j+1}"

def atom_values(engine, e, lam, dps):
    """values of all engine atoms; Li21 from lam cache (340-digit strings)."""
    mp.mp.dps = max(dps, 60)
    out = {}
    for s, key in engine.ATOMS.items():
        kind = key[0]
        if kind == 'log': out[s] = mp.log(mp.mpf(key[1]))
        elif kind == 'z2': out[s] = mp.pi**2/6
        elif kind == 'z3': out[s] = mp.zeta(3)
        elif kind == 'Li2': out[s] = mp.re(mp.polylog(2, mp.mpf(key[1])))
        elif kind == 'Li3': out[s] = mp.re(mp.polylog(3, mp.mpf(key[1])))
        elif kind == 'Li21':
            lbl = LAMKEY[(e, key[1], key[2])]
            out[s] = mp.mpf(lam[lbl])
        else: raise ValueError(key)
    return out

def Li3v(a):
    a = mp.mpf(a); return mp.re(mp.polylog(3, a))
def Li2v(a):
    a = mp.mpf(a); return mp.re(mp.polylog(2, a))

def rhs_terms_235():
    """edges 1-3: column semantics of sol vector (126)"""
    lists = pickle.load(open(O9B+'/Sunit_lists.pkl','rb'))
    Gargs = [sp.Rational(a) for a in lists['Gargs']]
    Dargs = [(sp.Rational(w), p) for w, p in lists['Dargs']]
    return Gargs, Dargs

def rhs_terms_2357():
    """edge 0: column semantics of sol vector (895)"""
    s0 = pickle.load(open(O9B+'/sol_B0.pkl','rb'))
    meta = s0['meta']
    sols = [sp.Rational(s) for s in meta['sols']]
    negw = [sp.Rational(s) for s in meta['negw']]
    logprods = [tuple(lp) for lp in meta['logprods']]
    return sols, negw, logprods

def eval_rhs_235(sol, Gargs, Dargs, dps):
    mp.mp.dps = dps
    Gv = [Li3v(g) for g in Gargs]
    Dv = [-Li2v(w)*mp.log(p) for w,p in Dargs]
    r = mp.fsum([mp.mpf(s_)*v for s_,v in zip(sol[:len(Gargs)], Gv)])
    r += mp.fsum([mp.mpf(s_)*v for s_,v in zip(sol[len(Gargs):], Dv)])
    return r

def eval_rhs_2357(sol, sols, negw, logprods, dps):
    mp.mp.dps = dps
    n = len(sols)
    r = mp.fsum([mp.mpf(s_)*Li3v(g) for s_,g in zip(sol[:n], sols)])
    Dv = []
    for w in negw:
        lw = Li2v(w)
        for p in (2,3,5,7):
            Dv.append(-lw*mp.log(p))
    r += mp.fsum([mp.mpf(s_)*v for s_,v in zip(sol[n:n+len(Dv)], Dv)])
    L = {p: mp.log(p) for p in (2,3,5,7)}
    r += mp.fsum([mp.mpf(s_)*L[a]*L[b]*L[c] for s_,(a,b,c) in zip(sol[n+len(Dv):], logprods)])
    return r
'''
open('oct9c_common.py','w').write(common)
print("written oct9c_common.py")
