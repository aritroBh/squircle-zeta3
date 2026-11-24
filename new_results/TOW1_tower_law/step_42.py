espace = '''

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
    return QSI(2*(d-2*m), 0, -2*d*m, 0) - QSI(4*(d-1)*(d-2*m),0,0,0).scale(Fr(1,2*d-3))
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
'''
open('towlib.py','a').write(espace)
import importlib, towlib
importlib.reload(towlib)
from towlib import *
print('e-space loaded')
print('e(2,1;0,0) =', e_off((2,1),0,0), 'vs rho', rho((2,1),0), rho_c((2,1),0))
print('w_slope(2,1) =', w_slope((2,1)))
print('data: u(2,1,k,2) slope should be sqrt2/4: from earlier: (k-3)*2/8 in p1: QSI(0,2/8,0,0) =', QSI(0,Fr(2,8)))
print('w_slope(3,0) =', w_slope((3,0)), ' (expect QSI(0,-2/24,0,0)= -sqrt2/12)')
print('ehat_22(0,0) =', ehat_pp(2,0,0), ' (expect stable 15sqrt2(1+i)/16)')
