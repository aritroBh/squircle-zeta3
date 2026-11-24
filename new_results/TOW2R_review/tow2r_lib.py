"""tow2r_lib.py -- TOW2R referee's independent machinery (from first principles).
Built from genintegrand.S_poly (master square root coefficient) + fastred.QSI exact arithmetic.
Nothing imported from tow2lib."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tow2_copy'))
from fractions import Fraction as Fr
import genintegrand as GI
from fastred import QSI, ZERO, ONE, IQ, AQ, BQ

def qsi_pmul(P, Q):
    out = [ZERO]*(len(P)+len(Q)-1)
    for a,ca in enumerate(P):
        if ca.is_zero(): continue
        for b,cb in enumerate(Q):
            if not cb.is_zero(): out[a+b] = out[a+b] + ca*cb
    while out and out[-1].is_zero(): out.pop()
    return out

def qsi_padd(P, Q):
    n = max(len(P), len(Q)); out = [ZERO]*n
    for k,c in enumerate(P): out[k]=out[k]+c
    for k,c in enumerate(Q): out[k]=out[k]+c
    while out and out[-1].is_zero(): out.pop()
    return out

def qsi_pscale(P, q): return [c*q for c in P]

def qsi_peval_shift(P, v0):
    """Q(y) = P(v0 + y)."""
    out = [ZERO]
    ypow = [ONE]
    for c in P:
        out = qsi_padd(out, qsi_pscale(ypow, c))
        ypow = qsi_pmul(ypow, [v0, ONE])
    return out

def Stilde(S):
    D = len(S)-1
    TV = [QSI(-1), QSI(2), QSI(1)]
    out = [ZERO]
    for s, c in enumerate(S):
        if c == 0: continue
        term = [ONE]
        for _ in range(s): term = qsi_pmul(term, TV)
        term = qsi_pscale(term, QSI(c))
        term = qsi_pscale(term, QSI(Fr(4)**(D-s)))
        term = [ZERO]*(D-s) + term
        out = qsi_padd(out, term)
    return out, D

def Phi_num_den(p, q):
    """Phi_p(v) = -(sqrt2/2) 8^d 4^{-D} v^{2d-1-D} Stilde (v^2+1)^{2-2d} / ((v^2-a^2)(v^2-b^2))
    as (numerator, denominator) QSI polynomials."""
    d = p+q
    S = GI.S_poly(d, p)
    St, D = Stilde(S)
    c = QSI(0, Fr(-(8**d), 2*4**D))     # -(sqrt2/2)*8^d/4^D
    num = qsi_pscale(St, c)
    den = [ZERO, ONE]                    # v (from v^{2d-1-D} = v^{-1})
    den = qsi_pmul(den, [ONE, ZERO, QSI(-6), ZERO, ONE])   # (v^2-a^2)(v^2-b^2) = v^4-6v^2+1
    v2p1 = [ONE, ZERO, ONE]
    for _ in range(2*d-2): den = qsi_pmul(den, v2p1)
    return num, den

def laurent_at_i(num, den, pole_order, nterms):
    """coeffs of y^{-pole_order}, ..., y^{-pole_order+nterms-1} of num(i+y)/den(i+y)."""
    Ny = qsi_peval_shift(num, IQ)
    Dy = qsi_peval_shift(den, IQ)
    k0 = 0
    while k0 < len(Dy) and Dy[k0].is_zero(): k0 += 1
    assert k0 == pole_order, (k0, pole_order)
    D1 = Dy[k0:]
    n = nterms
    d0inv = D1[0].inv()
    inv = [ZERO]*n
    inv[0] = d0inv
    for k in range(1, n):
        acc = ZERO
        for j in range(1, k+1):
            if j < len(D1): acc = acc + D1[j]*inv[k-j]
        inv[k] = -acc * d0inv
    R = [ZERO]*n
    for a in range(min(len(Ny), n)):
        if a >= len(Ny): break
        for b in range(n-a):
            R[a+b] = R[a+b] + Ny[a]*inv[b]
    return R   # R[t] = coeff of y^{-pole_order+t}

def rho_first_principles(p, q, delta_max):
    """u_{pq}^{(d)}(delta) for delta = 0..delta_max (initial tower, level k=d)."""
    d = p+q
    num, den = Phi_num_den(p, q)
    R = laurent_at_i(num, den, 2*d-2, delta_max+1)
    factor = Fr(2) if p > q else Fr(1)
    return [R[delta].scale(factor) for delta in range(delta_max+1)]

# ---------------- tau constants (TOW1 definitions) ----------------
def tau1(n):
    s = n+1
    return ((-1)**s) * ( (IQ-BQ)**(-(s+1)) + (IQ+AQ)**(-(s+1)) - IQ**(-(s+1)) )
def tau2(n):
    s = n+1
    return ((-1)**s) * ( (IQ-AQ)**(-(s+1)) + (IQ+BQ)**(-(s+1)) - IQ**(-(s+1)) )

# ---------------- T1 recursion (my implementation from TOW1's statement) ----------------
def level_towers(n_max, k, delta_max, rho_cache):
    """u_{pq}^{(k)}(delta) for all keys of degree <= n_max at level k, given
    initial data rho_cache[(p,q)] = [rho(0..delta_max)] for degrees n_max..k.
    Returns dict (p,q) -> [u(0..delta_max)]."""
    # towers_by_degree[d] = {(p,q): [u(0..delta_max)]}
    towers = {}
    towers[k] = {(p, k-p): rho_cache[(p, k-p)] for p in range(k, (k-1)//2, -1)}
    for d in range(min(k-1, n_max), 0, -1):
        pass
    # need descending d from k-1 down to 2; parents at degree d+1
    out = {}
    for d in range(k-1, 1, -1):
        if d > n_max: continue
        cur = {}
        for p in range(d, (d-1)//2, -1):
            q = d-p
            u = [ZERO]*(delta_max+1)
            # DL1 feed from (p+1,q)
            par = towers[d+1][(p+1, q)]
            for dp in range(0, delta_max+1):
                for dd_ in range(dp, delta_max+1):
                    pass
            for delta in range(0, delta_max+1):
                tot = ZERO
                for dp in range(0, delta+1):
                    c1 = tau1(delta-dp).scale(Fr(p+1, 2*d-dp-1))
                    tot = tot + par[dp]*c1
                    if p > q:
                        c2 = tau2(delta-dp).scale(Fr(q+1, 2*d-dp-1))
                        tot = tot + towers[d+1][(p, q+1)][dp]*c2
                if q == p-1:
                    tot = tot + fold_feed(p, towers[d+1][(p, p)], delta)
                u[delta] = tot
            cur[(p, q)] = u
        towers[d] = cur
        out.update(cur)
    return towers

from math import comb as _comb
def fold_feed(p, udiag, delta):
    """Fold_p(delta; u) = -p * sum_{mu=m}^{4p-4} (-i)^mu binom(mu-2, mu-m) i^{m-2} conj(g_mu),
    m = 4p-4-delta, g_mu = -sum_{dp} u(dp) tau1(4p-4-mu-dp)/(4p-3-dp)."""
    m = 4*p-4-delta
    tot = ZERO
    for mu in range(m, 4*p-3):
        g = ZERO
        for dp in range(0, 4*p-4-mu+1):
            if dp < len(udiag):
                g = g - udiag[dp]*tau1(4*p-4-mu-dp).scale(Fr(1, 4*p-3-dp))
        f = g.conj_i()
        tot = tot + f * (IQ**(-mu)) * QSI(_comb(mu-2, mu-m)) * (IQ**(m-2))
    return tot.scale(-p)
