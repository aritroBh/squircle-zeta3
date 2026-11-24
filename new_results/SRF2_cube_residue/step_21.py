import sys
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, ctx
from vlib import ts_nodes
flint.ctx.cap = 64
ctx.dps = 60; mp.mp.dps = 60

def tomp(x):
    return mp.mpf(x.str(70, radius=False, more=True))

def W_int(eps, h=arb(1)/32, T=4.6):
    """W(eps) = int int a^{2e-1} c^{e-1} [g - g0 - c*B1] da dc   (proper for e in (-1,0])."""
    nodes = ts_nodes(h, T)
    q = 1 + 2*eps
    tot = arb(0)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        g0 = (arb(-2*eps)*l1a).exp()/(1+a)*S0
        B1 = -q*a*(arb(-2*eps)*l1a).exp()/(1+a)**2*S0
        prefa = (arb(2*eps)*la).exp()/a
        row = arb(0)
        for (c, cm, wc, pc) in nodes:
            lc = c.log()
            C2 = c*c*(arb(-2*eps)*lc).exp()
            N = 1 + a + a*c
            g = (arb(-2*eps)*N.log()).exp()/N*(1+A2+A2*C2).sqrt()
            prefc = (arb(eps)*lc).exp()/c
            row += prefc*(g - g0 - c*B1)*wc
        tot += prefa*row*wa
    return tot

def Pc_int(eps, h=arb(1)/32, T=4.6):
    """Pc(eps) = -q/(e+1) * J(eps),  J = 1/(2e+1) + int_0^1 a^{2e}[(1+a)^{-1-q} S0 - 1] da."""
    nodes = ts_nodes(h, T)
    q = 1 + 2*eps
    J1 = arb(0)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(arb(-2*eps)*la).exp()
        S0 = (1+A2).sqrt()
        f = (arb(-2*eps)*l1a).exp()/(1+a)*S0 - 1
        J1 += (arb(2*eps)*la).exp()*f*wa
    J = 1/(2*eps+1) + J1
    return -q/(eps+1)*J, J

for e in [-0.9, -0.95, -0.99]:
    eps = arb(e)
    W = W_int(eps)
    Pc, J = Pc_int(eps)
    Vv = W + Pc
    print(f"eps={e}: W={mp.nstr(tomp(W),12)}  J={mp.nstr(tomp(J),12)}  V={mp.nstr(tomp(Vv),12)}  (1+eps)V={mp.nstr(tomp((1+eps)*Vv),15)}")
