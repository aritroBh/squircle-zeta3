"""cand.py -- candidate tower from the FORCED ascending recursion, evaluated numerically (exact Q(i)) at integer (d,m),
   base F(delta,0)(d,m) = a_delta(d,m) from MY first-principles engine (core.a_of); tau from MY DL1/DL2 series."""
import sys
sys.dont_write_bytecode=True
from core import *
from math import comb

_F = {}
def reset():
    _F.clear()
def Fnum(delta, s, D, M, tau1=None, tau2=None, base=None):
    """F(delta,s)(D,M) in Q(i).  tau1/tau2/base can be overridden (mutation tests)."""
    tau1 = tau1 or TAU1; tau2 = tau2 or TAU2
    if s > delta // 2: return ZERO
    key = (delta, s, D, M, id(tau1), id(tau2), id(base))
    if key in _F: return _F[key]
    if s == 0:
        r = (base(D, M)[delta] if base else a_of(D, M)[delta])
    else:
        r = Lnum(delta, s - 1, D, M, tau1, tau2, base) - Fnum(delta, s - 1, D, M, tau1, tau2, base)
    _F[key] = r
    return r

def Lnum(delta, s, D, M, tau1, tau2, base):
    """Lambda^{(delta)}(s)(D,M) = i(2D-1) sum_{dp<=delta} [tau1_{delta-dp} F(dp,s)(D+1,M+1) - tau2_{delta-dp} F(dp,s)(D+1,M)]/(2D-dp-1)"""
    tot = ZERO
    for dp in range(delta + 1):
        fA = Fnum(dp, s, D + 1, M + 1, tau1, tau2, base)
        fB = Fnum(dp, s, D + 1, M, tau1, tau2, base)
        tot = tot + (tau1[delta - dp] * fA - tau2[delta - dp] * fB) * fmpq(1, 2 * D - dp - 1)
    return tot * (I * (2 * D - 1))

def closing(delta, D, M, **kw):
    """F(delta, floor(delta/2)+1)(D,M) = Lam(delta,top) - F(delta,top)"""
    top = delta // 2
    tau1 = kw.get('tau1') or TAU1; tau2 = kw.get('tau2') or TAU2; base = kw.get('base')
    return Lnum(delta, top, D, M, tau1, tau2, base) - Fnum(delta, top, D, M, tau1, tau2, base)

def cand_value(p, q, delta, k, tau1=None, tau2=None, base=None, diag_rule='stable'):
    """candidate u~(delta) for key (p,q) at level k (k>=d off-diag, k>=2p+1 diag)."""
    tau1 = tau1 or TAU1; tau2 = tau2 or TAU2
    d = p + q
    r0 = rho0_tilde(p, q)
    if p > q:
        tot = ZERO
        for s in range(delta // 2 + 1):
            tot = tot + Fnum(delta, s, d, p, tau1, tau2, base) * comb(k - d, s)
        return r0 * tot
    else:
        # diagonal (p,p): ehat(delta,s) = (p+1) sum_{dp<=delta} tau1_{delta-dp} e_par(dp,s)/(4p-dp-1), e_par = rho~0(p+1,p) F(dp,s)(2p+1,p+1)
        rp = rho0_tilde(p + 1, p)
        tot = ZERO
        for s in range(delta // 2 + 1):
            eh = ZERO
            for dp in range(delta + 1):
                eh = eh + tau1[delta - dp] * Fnum(dp, s, 2 * p + 1, p + 1, tau1, tau2, base) * fmpq(1, 4 * p - dp - 1)
            eh = eh * rp * (p + 1)
            n = k - 2 * p - 1
            tot = tot + eh * comb(n, s)
        return tot
