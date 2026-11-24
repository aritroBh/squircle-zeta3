lib = '''
import sys
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, ctx
from vlib import ts_nodes

def tomp(x):
    return mp.mpf(x.str(120, radius=False, more=True))

def expm1m(z):
    """expm1(z) - z = sum_{m>=2} z^m/m! (stable for small z)."""
    if abs(z) > arb('0.5'):
        return z.exp() - 1 - z
    t = z*z/2
    s = t
    for m in range(3, 200):
        t = t*z/m
        s = s + t
        if abs(t) < arb(10)**(-(ctx.dps+15)):
            break
    return s

def u_minus_log1p(u):
    """u - log(1+u), stable for small u."""
    if u > arb('0.25'):
        return u - (1+u).log()
    t = u*u/2
    s = t
    sgn = -1
    uj = u*u
    for j in range(3, 600):
        uj = uj*u
        t = uj/j
        s = s + sgn*t
        sgn = -sgn
        if abs(t) < arb(10)**(-(ctx.dps+15)):
            break
    return s

def Phi(u, q):
    """(1+u)^{-q} - 1 + q*u, stable."""
    if u > arb('0.25'):
        return (arb(-q)*(1+u).log()).exp() - 1 + q*u
    z = -q*(1+u).log()
    return expm1m(z) + q*u_minus_log1p(u)

def W_int(eps, h=arb(1)/32, T=4.6):
    """W(eps) = int int_[0,1]^2 a^{2e-1} c^{e-1} [g - g0 - c*B1] da dc  (proper for e in (-1,0]).
       bracket = (1+a)^{-1-2e} [ Phi(u) S + (1 - q u) * (A2*C2)/(S+S0) ],  q=1+2e, u=ac/(1+a)."""
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
        outfact = (arb(-2*eps)*l1a).exp()/(1+a)   # (1+a)^{-1-2e}
        tot_row = arb(0)
        for (c, C2, prefc, wc) in cpre:
            u = a*c/(1+a)
            S = (1+A2+A2*C2).sqrt()
            br = Phi(u, q)*S + (1 - q*u)*(A2*C2)/(S+S0)
            tot_row += prefc*br*wc
        tot += prefa*outfact*tot_row*wa
    return tot

def J_int(eps, h=arb(1)/32, T=4.6):
    """J(eps) = int_0^1 a^{2e} (1+a)^{-2-2e} S0 da = 1/(2e+1) + J1 (J1 proper for e>-1)."""
    nodes = ts_nodes(h, T)
    J1 = arb(0)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        f = S0*expm1m(arb(-(2+2*eps))*l1a) + A2/(S0+1)
        J1 += (arb(2*eps)*la).exp()*f*wa
    return 1/(2*eps+1) + J1

def Pc_int(eps, h=arb(1)/32, T=4.6):
    q = 1 + 2*eps
    return -q/(eps+1)*J_int(eps, h, T)

def V_cont(eps, h=arb(1)/32, T=4.6):
    """V(eps) for eps in (-1, 0]:  V = W + Pc."""
    return W_int(eps, h, T) + Pc_int(eps, h, T)

def Jstar_int(eps, h=arb(1)/32, T=4.6):
    """J*(eps) = int_0^1 a^{2e} S0 da = 1/(2e+1) + int a^{2e}(S0-1) da (planner's regularised)."""
    nodes = ts_nodes(h, T)
    J1 = arb(0)
    for (a, am, wa, pa) in nodes:
        la = a.log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        J1 += (arb(2*eps)*la).exp()*(A2/(S0+1))*wa
    return 1/(2*eps+1) + J1

def I_quad(h=arb(1)/32, T=4.6):
    """I = int_0^1 (sqrt(1+a^4)-1)/a^2 da  = int_0^1 a^2/( (sqrt(1+a^4)+1) ) da... stable."""
    nodes = ts_nodes(h, T)
    s = arb(0)
    for (a, am, wa, pa) in nodes:
        S0 = (1+a**4).sqrt()
        s += a*a/(S0+1)*wa   # (sqrt(1+a^4)-1)/a^2 = a^2/(S0+1)
    return s
'''
with open('<project>/research/08_theory/SRF2_cube_residue/srf2_lib.py', 'w') as f:
    f.write(lib)
print("written", len(lib))
