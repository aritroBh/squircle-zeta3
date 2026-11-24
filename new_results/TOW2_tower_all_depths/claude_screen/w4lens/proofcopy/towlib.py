
"""towlib: tower data at v=+i, depth-indexed, for the TOW1 job."""
import sys, os, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
WTS = os.path.normpath(os.path.join(HERE, '..', 'WTS_weight_k_structure', 'theory'))
from fractions import Fraction as Fr
from fastred import *

def q4(q): return (q.p0, q.p1, q.p2, q.p3)
def unq4(t): return QSI(*t)

_towers = {}
def tower(k):
    if k not in _towers:
        _towers[k] = pickle.load(open(os.path.join(WTS, f'tower_{k}.pkl'), 'rb'))
    return _towers[k]

def hcoef(i, j, k, r):
    """coefficient of (v-i)^{-r} in HigherPart(W_{ij}^{(k)}) at +i, as QSI; 0 if absent."""
    H = tower(k).get((i,j))
    if H is None: return ZERO
    return unq4(H.get(r, (0,0,0,0)))

def deg(key): return key[0]+key[1]

def u(i, j, k, delta):
    """depth-delta tower coefficient: order r = 2d-2-delta."""
    d = i+j
    r = 2*d-2-delta
    if r < 2: raise ValueError('delta too large')
    return hcoef(i, j, k, r)

_binteg = {}
def rho(key, delta):
    """initial data: depth-delta coefficient of the level-d integrand (HigherPart of R_ij)."""
    d = deg(key)
    if d not in _binteg:
        _binteg[d] = pickle.load(open(os.path.join(WTS, f'binteg_{d}.pkl'), 'rb'))
    P = _binteg[d][key]
    # P is a PF [real, ipos, ineg, const] of QSI (as loaded); higher part at +i: ipos[1:]
    ipos = P[1]
    r = 2*d-2-delta
    if r-1 < len(ipos):
        return ipos[r-1]
    return ZERO

def keys_of_degree(d):
    return [(i, d-i) for i in range(d, (d-1)//2, -1)]

# Taylor coefficients of DL1, DL2 at v=i: tau1[n] = DL1^{(n+1)}(i)/(n+1)!
# DL1 = 1/(v-b)+1/(v+a)-1/v  =>  DL1^{(s)}(i)/s! = (-1)^s[1/(i-b)^{s+1}+1/(i+a)^{s+1}-1/i^{s+1}]
def tau1(n):
    s = n+1   # tau1[n] = DL1^{(s)}(i)/s! = (-1)^s sum_r (i-r)^{-(s+1)}
    return ((-1)**s) * ( (IQ-BQ)**(-(s+1)) + (IQ+AQ)**(-(s+1)) - IQ**(-(s+1)) )
def tau2(n):
    s = n+1
    return ((-1)**s) * ( (IQ-AQ)**(-(s+1)) + (IQ+BQ)**(-(s+1)) - IQ**(-(s+1)) )


def fold_feed(i, k, delta):
    """the +i order-(2d-2-delta) coefficient of -i*F(Gm_ii DL1) at level k (d=2i-1), from tower data.
       = -i * sum_{mu=m}^{4i-4} f_mu (-i)^mu binom(mu-2, mu-m) i^{m-2}, f_mu = conj(g_mu),
       g_mu = -sum_{dp=0}^{4i-4-mu} u_ii(k,dp) tau1[4i-4-mu-dp]/(4i-3-dp)."""
    d = 2*i-1
    m = 4*i-4-delta
    tot = ZERO
    from math import comb
    for mu in range(m, 4*i-3):   # mu <= 4i-4
        g = ZERO
        for dp in range(0, 4*i-4-mu+1):
            g = g - u(i,i,k,dp) * tau1(4*i-4-mu-dp).scale(Fr(1, 4*i-3-dp))
        f = g.conj_i()
        tot = tot + f * (IQ**(-mu)) * Fr(comb(mu-2, mu-m)) * (IQ**(m-2))
    return tot.scale(-i)   # W-term: -(integer i) * phi_m

def Lambda_ij(i, j, k, delta, U):
    """the recursion feed (k>d part): value of (Lambda_d U_{d+1})(i,j) at depth delta,
       U(key, dp) = degree-(d+1) depth-dp data at level k."""
    d = i+j
    tot = ZERO
    for dp in range(0, delta+1):
        c1 = tau1(delta-dp).scale(Fr(i+1, 2*d-dp-1))
        tot = tot + U((i+1,j), dp) * c1
        if i > j:
            c2 = tau2(delta-dp).scale(Fr(j+1, 2*d-dp-1))
            tot = tot + U((i,j+1), dp) * c2
    if j == i-1:
        # fold term: need g_mu built from U((i,i), dp)
        from math import comb
        m = 4*i-4-delta
        ft = ZERO
        for mu in range(m, 4*i-3):
            g = ZERO
            for dp in range(0, 4*i-4-mu+1):
                g = g - U((i,i), dp) * tau1(4*i-4-mu-dp).scale(Fr(1, 4*i-3-dp))
            ft = ft + g.conj_i() * (IQ**(-mu)) * Fr(comb(mu-2, mu-m)) * (IQ**(m-2))
        tot = tot + ft.scale(-i)   # -(integer i) * phi_m
    return tot


def polyval(key, delta, k, diag_shift=None):
    """tower-law polynomial for key at depth delta, evaluated at level k.
       Fit from levels k0..k0+deg where k0 = d (off-diag) or d+1 (diagonal)."""
    d = key[0]+key[1]
    if diag_shift is None:
        diag_shift = (key[0]==key[1])
    k0 = d + (1 if diag_shift else 0)
    degf = delta//2
    xs = list(range(k0, k0+degf+1))
    ys = [u(key[0],key[1],x,delta) for x in xs]
    out = ZERO
    for a,xa in enumerate(xs):
        La = Fr(1)
        for b,xb in enumerate(xs):
            if b!=a: La *= Fr(k-xb, xa-xb)
        out = out + ys[a].scale(La)
    return out

def feed_diag_into_child(p, k, delta, Udiag):
    i = p; j = p-1
    return Lambda_ij(i,j,k,delta, lambda kk,dp: Udiag(dp) if kk==(p,p) else ZERO)




# ============ Newton-coefficient families e(delta,s) for the tower law, depths <= 2 ============
from genintegrand import S_poly, binom_half
from math import comb
def C(n,k): return comb(n,k) if (0 <= k <= n) else 0
_tstar = QSI(Fr(1,2),0,Fr(1,2),0)
def S_tstar(d,m):
    return QSI(Fr(binom_half(d))) * (IQ**d) * Fr((-1)**m * C(d,m))
def Sp_tstar(d,m):
    Bd_prime = QSI(Fr(C(d,m))) * (IQ**d) * Fr((-1)**m) * QSI(2*(d-2*m),0,-2*d,0)
    c = Fr(C(d-1,m)*(d-m-1) - C(d-1,m-1)*(m-1))
    Bdm1 = QSI(c) * (IQ**(d-1)) * Fr((-1)**m)
    return QSI(Fr(binom_half(d))) * Bd_prime + QSI(Fr(binom_half(d-1))) * (2*IQ) * Bdm1
def r_i(d, m):
    """(S_m''/S_m)(t*) = 2(d-2m) - 2 d m_imag? NO: = 2(d-2m) - 2*d*m*IUNIT ... see formula."""
    return QSI(2*(d-2*m), 0, -2*d, 0) - QSI(4*(d-1)*(d-2*m),0,0,0).scale(Fr(1,2*d-3))
    # r_m = 2(d-2m) - 2 d m * i - 4(d-1)(d-2m)/(2d-3)

def degS(d, m):
    P = S_poly(d,m)
    return len(P)-1 if P else 0

def A_i_value(d, i):
    D = degS(d, i)
    St_i = (QSI(4)*IQ)**D * S_tstar(d, i)
    denom = (IQ**2 - AQ**2)*(IQ**2 - BQ**2)
    C0 = QSI(0, Fr(-1,2)) * QSI(8**d) / QSI(4**D)
    return C0 * (IQ**(2*d-1-D)) * St_i * ((2*IQ)**(2-2*d)) / denom

def rho_c(key, delta):
    """closed-form initial data, depths 0,1,2 (single piece times (1+[i>j]))."""
    i, j = key; d = i+j
    mult = 2 if i > j else 1
    A = A_i_value(d, i)
    if delta == 0: return A.scale(mult)
    if delta == 1: return (A * QSI(0,0,-(d-2),0)).scale(mult)
    if delta == 2:
        rr = r_i(d, i)
        val = QSI(Fr(3*d-3,2)) - rr*(IQ.scale(Fr(1,2))) - QSI((d-2)**2)
        return (A*val).scale(Fr(mult,2))

# e-coefficients (off-diagonal): e(0,0),e(1,0),e(2,0) = rho; e(2,1) = slope w (computed below)
def e_off(key, delta, s):
    if (delta, s) == (0,0): return rho_c(key, 0)
    if (delta, s) == (1,0): return rho_c(key, 1)
    if (delta, s) == (2,0): return rho_c(key, 2)
    if (delta, s) == (2,1): return w_slope(key)
    return ZERO

def ehat_pp(p, delta, s):
    """diagonal stable Newton coefficients from off-diagonal parent (p+1,p)."""
    if delta == 0 and s == 0:
        return e_off((p+1,p), 0, 0) * tau1(0).scale(Fr(p+1, 4*p-1))
    if delta == 1 and s == 0:
        return ( e_off((p+1,p),1,0) * tau1(0).scale(Fr(p+1, 4*p-2))
               + e_off((p+1,p),0,0) * tau1(1).scale(Fr(p+1, 4*p-1)) )
    if delta == 2 and s == 0:
        return ( e_off((p+1,p),2,0) * tau1(0).scale(Fr(p+1, 4*p-3))
               + e_off((p+1,p),1,0) * tau1(1).scale(Fr(p+1, 4*p-2))
               + e_off((p+1,p),0,0) * tau1(2).scale(Fr(p+1, 4*p-1)) )
    if delta == 2 and s == 1:
        return w_slope((p+1,p)) * tau1(0).scale(Fr(p+1, 4*p-3))
    return ZERO

def e_any(key, delta, s):
    i, j = key
    if i == j: return ehat_pp(i, delta, s)
    return e_off(key, delta, s)

def w_slope(key):
    """e(2,1) = [depth-2 feed of e(.,0)] - rho(2). Computed via Lambda on a fake level."""
    i, j = key
    feed = LambdaE_ij(i, j, 2, 0)   # depth 2, Newton level sigma=0
    return feed - rho_c(key, 2)

def LambdaE_ij(i, j, delta, sigma):
    """Newton-space recursion RHS at (delta, sigma), WITHOUT LHS terms:
       sum_{delta'<=delta} [ (i+1) tau1_{delta-delta'}/(2d-delta'-1) Ep1 + [i>j](j+1) tau2 ... Ep2 ] + fold.
       EpX(delta',sigma): off-diag parent: e(parent, delta', sigma); diagonal parent: alt-sum_t (-1)^t e_hat(delta', sigma+t)."""
    d = i+j
    tot = ZERO
    for dp in range(0, delta+1):
        c1 = tau1(delta-dp).scale(Fr(i+1, 2*d-dp-1))
        tot = tot + e_any((i+1,j), dp, sigma) * c1
        if i > j:
            c2 = tau2(delta-dp).scale(Fr(j+1, 2*d-dp-1))
            tot = tot + Ediag_conv((i,j+1), dp, sigma) * c2
    if j == i-1:
        # fold: -i_int * phi, phi from conj(g_mu), g_mu from diagonal u_ii via tau1;
        # in Newton space: diagonal entries pass through alt-sum conversion
        from math import comb as _comb
        m = 4*i-4-delta
        ft = ZERO
        for mu in range(m, 4*i-3):
            g = ZERO
            for dp in range(0, 4*i-4-mu+1):
                g = g - Ediag_conv((i,i), dp, sigma) * tau1(4*i-4-mu-dp).scale(Fr(1, 4*i-3-dp))
            ft = ft + g.conj_i() * (IQ**(-mu)) * Fr(_comb(mu-2, mu-m)) * (IQ**(m-2))
        tot = tot + ft.scale(-i)
    return tot

def Ediag_conv(key, delta, sigma):
    """if off-diagonal: e; if diagonal (p,p): alt sum_{t>=0} (-1)^t ehat_pp(delta, sigma+t)."""
    i, j = key
    if i != j: return e_off(key, delta, sigma)
    p = i
    out = ZERO
    t = 0
    while True:
        e = ehat_pp(p, delta, sigma+t)
        # ehat is zero for sigma+t > floor(delta/2); stop when all further terms vanish
        if sigma+t > delta//2: break
        out = out + e.scale((-1)**t)
        t += 1
    return out


def e_off(key, delta, s):
    if (delta, s) == (0,0): return rho_c(key, 0)
    if (delta, s) == (1,0): return rho_c(key, 1)
    if (delta, s) == (2,0): return rho_c(key, 2)
    if (delta, s) == (2,1):
        i,j = key; d = i+j
        return rho_c(key, 0).scale(Fr(-1, 2*d-3))     # closed-form slope -rho(0)/(2d-3)
    return ZERO

def ehat_pp(p, delta, s):
    if delta == 0 and s == 0:
        return e_off((p+1,p), 0, 0) * tau1(0).scale(Fr(p+1, 4*p-1))
    if delta == 1 and s == 0:
        return ( e_off((p+1,p),1,0) * tau1(0).scale(Fr(p+1, 4*p-2))
               + e_off((p+1,p),0,0) * tau1(1).scale(Fr(p+1, 4*p-1)) )
    if delta == 2 and s == 0:
        return ( e_off((p+1,p),2,0) * tau1(0).scale(Fr(p+1, 4*p-3))
               + e_off((p+1,p),1,0) * tau1(1).scale(Fr(p+1, 4*p-2))
               + e_off((p+1,p),0,0) * tau1(2).scale(Fr(p+1, 4*p-1)) )
    if delta == 2 and s == 1:
        return e_off((p+1,p),2,1) * tau1(0).scale(Fr(p+1, 4*p-3))
    return ZERO


def U_e(key, k, delta):
    """tower-law candidate from e-coefficients; diagonal at k=2p uses rho (initial)."""
    i, j = key; d = i+j
    if delta > 2*d-4: return None
    if i == j:
        p = i
        if k == 2*p:
            return rho_c(key, delta)
        return sum( (ehat_pp(p, delta, s).scale(comb(k-2*p-1, s)) for s in range(0, delta//2+1)), ZERO )
    return sum( (e_off(key, delta, s).scale(comb(k-d, s)) for s in range(0, delta//2+1)), ZERO )

def check_tower_law(kmax=14, dmax=10):
    bad = []
    for k in range(2, kmax+1):
        for d in range(2, min(dmax, k)+1):
            for key in keys_of_degree(d):
                i, j = key
                for delta in range(0, min(2, 2*d-4)+1):
                    ue = U_e(key, k, delta)
                    actual = u(i,j,k,delta)
                    onlaw = (i > j) or (k >= 2*i+1)
                    if onlaw and ue != actual:
                        bad.append((key,k,delta,ue,actual))
    return bad

def check_newton(dmax=9):
    """k-free Newton identities: e(delta,sigma)+e(delta,sigma+1) == LambdaE_ij(i,j,delta,sigma),
       for off-diagonal keys; diagonal keys: ehat == parent formula (by construction, but check vs data)."""
    bad = []
    for d in range(2, dmax+1):
        for key in keys_of_degree(d):
            i, j = key
            if i == j: continue
            for delta in (0,1,2):
                if delta > 2*d-4: continue
                for sigma in range(0, delta//2+1):
                    lhs = e_off(key, delta, sigma) + e_off(key, delta, sigma+1)
                    rhs = LambdaE_ij(i, j, delta, sigma)
                    if lhs != rhs:
                        bad.append((key,delta,sigma,lhs,rhs))
    return bad
