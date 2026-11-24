import sys, os, time
work = '<project>/research/08_theory/TOW2_tower_all_depths/claude_screen/scripts_rerun'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
from towlib import *
from fractions import Fraction as Fr
from math import comb

# closed-form e-family as QSI values: e(delta,s) for key, via tow2lib formulas
d_sym, m_sym = T.d, T.m
def qsi_of(expr, dd, mm):
    v = sp.expand(expr.subs({d_sym: dd, m_sym: mm}))
    re = sp.simplify(sp.re(v)); im = sp.simplify(sp.im(v))
    assert not re.has(sp.sqrt(2)) and not im.has(sp.sqrt(2))
    return QSI(Fr(sp.Rational(re)), 0, Fr(sp.Rational(im)), 0)

def e_off(key, delta, s):
    """e(delta,s) absolute = closed form * e00(key)."""
    i_, j_ = key; dd = i_+j_
    if s > delta//2: return ZERO
    if s == 0: return rho_c(key, delta) if delta <= 2 else rho(key, delta)
    # s >= 1: use closed forms
    if delta == 2 and s == 1: form = -1/(2*dd_sym-3)
    elif delta == 3 and s == 1: form = sp.I*(dd_sym-3)/(2*dd_sym-3)
    elif delta == 4 and s == 1:
        form = (dd_sym**3 - 9*dd_sym**2 + sp.Rational(113,4)*dd_sym - sp.I*dd_sym/2 + sp.I*mm_sym - 29)/((2*dd_sym-3)*(2*dd_sym-5))
    elif delta == 4 and s == 2: form = sp.Rational(2)/((2*dd_sym-3)*(2*dd_sym-5))
    else: raise ValueError((delta,s))
    return rho(key, 0) * qsi_of(form, dd, i_)

def ehat_pp(p_, delta, s):
    """ehat_pp(delta,s) absolute: (p+1) sum_dp' tau1 e_{p+1,p}(dp',s)/(4p-dp'-1)."""
    tot = ZERO
    for dpp in range(0, delta+1):
        tot = tot + e_off((p_+1,p_), dpp, s) * tau1(delta-dpp).scale(Fr(p_+1, 4*p_-dpp-1))
    return tot

def U_e(key, k, delta):
    i_, j_ = key; dd = i_+j_
    if delta > 2*dd-4: return None
    if i_ == j_:
        if k == 2*i_: return rho(key, delta)
        return sum( (ehat_pp(i_, delta, s).scale(comb(k-2*i_-1, s)) for s in range(0, delta//2+1)), ZERO )
    return sum( (e_off(key, delta, s).scale(comb(k-dd, s)) for s in range(0, delta//2+1)), ZERO )

bad = []; tested = 0
t0 = time.time()
for k in range(2, 15):
    for dd in range(2, min(11, k)+1):
        for key in keys_of_degree(dd):
            i_, j_ = key
            for delta in range(0, min(4, 2*dd-4)+1):
                ue = U_e(key, k, delta)
                actual = u(i_, j_, k, delta)
                onlaw = (i_ > j_) or (k >= 2*i_+1)
                tested += 1
                if onlaw and ue != actual:
                    bad.append((key,k,delta))
print(f'tower law delta<=4 vs data: {tested} cells, {len(bad)} bad  ({round(time.time()-t0,1)}s)')
if bad[:5]: print(bad[:5])
