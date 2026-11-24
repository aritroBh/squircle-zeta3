import sys
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, ctx
from vlib import ts_nodes
flint.ctx.cap = 64
ctx.dps = 80; mp.mp.dps = 80

def tomp(x):
    return mp.mpf(x.str(90, radius=False, more=True))

def expm1m(z):
    """expm1(z) - z = sum_{m>=2} z^m/m! (stable)."""
    if abs(z) > arb('0.5'):
        return z.exp() - 1 - z
    t = z*z/2
    s = t
    for m in range(3, 80):
        t = t*z/m
        s = s + t
        if abs(t) < arb(10)**(-95):
            break
    return s

def u_minus_log1p(u):
    """u - log(1+u), stable for small u."""
    if u > arb('0.25'):
        return u - (1+u).log()
    # series sum_{j>=2} (-1)^j u^j / j
    t = u*u/2
    s = t
    sgn = -1
    uj = u*u
    for j in range(3, 400):
        uj = uj*u
        t = uj/j
        s = s + sgn*t
        sgn = -sgn
        if abs(t) < arb(10)**(-95):
            break
    return s

def Phi(u, q):
    """(1+u)^{-q} - 1 + q*u, stable."""
    if u > arb('0.25'):
        return (arb(-q)*(1+u).log()).exp() - 1 + q*u
    z = -q*(1+u).log()
    return expm1m(z) + q*u_minus_log1p(u)

def W_int(eps, h=arb(1)/32, T=4.6):
    """W(eps) = int int a^{2e-1} c^{e-1} [g - g0 - c*B1] da dc, proper for e in (-1,0].
       bracket = (1+a)^{-1-q} [ Phi(u) S + (1 - q u) * (A2*C2)/(S+S0) ]."""
    q = 1 + 2*eps
    nodes = ts_nodes(h, T)
    cpre = []
    for (c, cm, wc, pc) in nodes:
        lc = c.log()
        C2 = c*c*(arb(-2*eps)*lc).exp()
        prefc = (arb(eps)*lc).exp()/c
        cpre.append((c, C2, prefc, wc))
    tot = arb(0)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        prefa = (arb(2*eps)*la).exp()/a
        outfact = (arb(-2*eps)*l1a).exp()/(1+a)   # (1+a)^{-1-2eps} = (1+a)^{-q}? NO:
        tot_row = arb(0)
        for (c, C2, prefc, wc) in cpre:
            u = a*c/(1+a)
            S = (1+A2+A2*C2).sqrt()
            br = Phi(u, q)*S + (1 - q*u)*(A2*C2)/(S+S0)
            tot_row += prefc*br*wc
        tot += prefa*outfact*tot_row*wa
    return tot

def J_int(eps, h=arb(1)/32, T=4.6):
    """J(eps) = int_0^1 a^{2e} (1+a)^{-2-2e} S0 da  = 1/(2e+1) + J1 (proper J1 for e>-1)."""
    nodes = ts_nodes(h, T)
    J1 = arb(0)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        f = S0*expm1m(arb(-(2+2*eps))*l1a) + A2/(S0+1)
        J1 += (arb(2*eps)*la).exp()*f*wa
    return 1/(2*eps+1) + J1

for e in [-0.9, -0.95, -0.99]:
    eps = arb(e)
    W = W_int(eps)
    J = J_int(eps)
    q = 1 + 2*eps
    Pc = -q/(eps+1)*J
    Vv = W + Pc
    print(f"eps={e}: W={mp.nstr(tomp(W),12)}  J={mp.nstr(tomp(J),12)}  V={mp.nstr(tomp(Vv),12)}  (1+eps)V={mp.nstr(tomp((1+eps)*Vv),15)}")
