"""candlib.py -- candidate tower built ONLY from closed-form e-family (sympy exprs) + tau table; T1 pointwise test driver."""
import sys, time
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from math import comb
import sympy as sp
import tow2lib as T2, tow3lib as T3
from t1lib import *

def gq(z):
    z = sp.nsimplify(z); return GQ(Fr(int(sp.re(z).p), int(sp.re(z).q)), Fr(int(sp.im(z).p), int(sp.im(z).q)))
def make_taus(tab=None):
    t1 = [gq(T2.tau1f(n)) for n in range(8)]; t2 = [gq(T2.tau2f(n)) for n in range(8)]
    return Taus(t1, t2)
TAUS = make_taus()

def base_exprs():
    ex = {}
    for delta in range(0, 7):
        for s in range(0, delta//2+1):
            ex[(delta, s)] = T3.e(delta, s)
    return ex

class Cand:
    def __init__(self, exprs, taus=TAUS, diag_initial_shift=None, ehat_mod=None):
        self.fam = Family(exprs); self.taus = taus
        self._eh = {}
        self.diag_initial_shift = diag_initial_shift   # mutation: callable(p,delta)->QSI added to rho_pp
        self.ehat_mod = ehat_mod
    def E(self, p, q, delta, s):               # absolute e_{pq}(delta,s) for an off-diagonal key (QSI)
        return E00(p, q).mul(self.fam(delta, s, p+q, p))
    def Ehat(self, p, delta, s):               # absolute ehat_pp(delta,s)
        key = (p, delta, s)
        r = self._eh.get(key)
        if r is None:
            tot = ZERO
            for dpp in range(delta+1):
                c = self.taus.t1[delta-dpp] * Fr(p+1, 4*p-dpp-1)
                tot = tot + self.E(p+1, p, dpp, s).mul(c)
            r = tot; self._eh[key] = r
        return r
    def get(self, p, q, k, delta):
        d = p+q
        assert delta <= 2*d-4 and delta <= 6, (p, q, k, delta)
        if p > q:
            assert k >= d
            tot = ZERO
            for s in range(delta//2+1):
                tot = tot + self.E(p, q, delta, s).mul(comb(k-d, s) if True else 0) if False else tot + self.E(p, q, delta, s).mul(GQ(comb(k-d, s)))
            return tot
        # diagonal
        assert k >= 2*p
        if k == 2*p:
            r = E00(p, p).mul(self.fam(delta, 0, 2*p, p))
            if self.diag_initial_shift is not None: r = r + self.diag_initial_shift(p, delta)
            return r
        tot = ZERO
        for s in range(delta//2+1):
            tot = tot + self.Ehat(p, delta, s).mul(GQ(comb(k-2*p-1, s)))
        return tot

def t1_test(cand, dlist, kspan=6, deltas=range(0, 7), stop_after=None, fold_sign=1, binom_fn=None, taus_for_T1=None):
    """test candidate against T1 at keys of degree in dlist, k=d+1..d+kspan.  returns (tested, bad list, per-delta count)"""
    tested = 0; bad = []; per = {}
    taus = taus_for_T1 or cand.taus
    for dd in dlist:
        for p in range(dd, (dd-1)//2, -1):   # p >= q
            q = dd-p
            for k in range(dd+1, dd+kspan+1):
                if p == q and k == 2*p: pass
                for delta in deltas:
                    if delta > 2*dd-4: continue
                    if p == q and k == 2*p: continue            # diagonal initial level is data, not T1
                    lhs = cand.get(p, q, k, delta)
                    rhs = T1rhs(cand.get, taus, p, q, k, delta, fold_sign=fold_sign, binom_fn=binom_fn)
                    tested += 1; per[delta] = per.get(delta, 0)+1
                    if not (lhs == rhs):
                        bad.append((p, q, k, delta))
                        if stop_after and len(bad) >= stop_after: return tested, bad, per
    return tested, bad, per
