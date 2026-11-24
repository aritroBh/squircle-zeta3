
import sys, time
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, arb_series, ctx
from vlib import ts_nodes, ser, tompf, V_series
flint.ctx.cap = 32
ctx.dps = 50; mp.mp.dps = 50

def Vcont_series(eps0, nlen=2, h=arb(1)/32, T=4.2):
    X = ser(eps0, nlen)
    nodes = ts_nodes(h, T)
    cpre = [(c, c*c*(X*(-2*c.log())).exp(), (X*c.log()).exp()/c, wc) for (c,cm,wc,pc) in nodes]
    full = arb_series([0], prec=nlen)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(X*(-2*la)).exp()
        S0 = (1+A2).sqrt()
        g0 = (X*(-2*l1a)).exp()/(1+a)*S0
        B1 = -(1+2*X)*a*(X*(-2*l1a)).exp()/(1+a)**2*S0
        prefa = (X*(2*la)).exp()/a
        pa_factor = (X*(-(1+2*X)*0))  # dummy
        row = arb_series([0], prec=nlen)
        for (c, C2, prefc, wc) in cpre:
            u = a*c/(1+a)
            S = (1+A2+A2*C2).sqrt()
            # stable: g - g0 = (1+a)^{-1-2e}[ ((1+u)^{-1-2e}-1) S + (S-S0) ]
            T1 = (X*(-(1)*u.log1p())).exp()  # (1+u)^{-e} ... need -(1+2e)*log1p(u)
            T1 = (X*(-1)*u.log1p() - 2*X*u.log1p()).exp() if False else None
            # simpler: ((1+u)^{-1-2e} - 1) = expm1((-1-2e)*log1p(u))
            arg = (-1-2*X)*u.log1p()
            T1 = arg.expm1()*S
            T2 = (A2*C2)/(S+S0)
            gg = (X*(-2*l1a)).exp()/(1+a)*(T1 + T2)
            A1 = -(1+2*X)*c
            row += prefc*(gg - a*A1 - c*B1)*wc
        full += prefa*row*wa
    Pa = -(1+2*X)/((2*X+1)*(X+1))
    Pc = arb_series([0], prec=nlen)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(X*(-2*la)).exp()
        Ic = (X*(2*la)).exp()*(X*(-2*l1a)).exp()/(1+a)**2*(1+A2).sqrt()
        Pc += -(1+2*X)/(X+1)*Ic*wa
    return full + Pa + Pc

# sanity at eps=+0.1 vs raw V_series (raw returns list of coeffs at 0)
vf, _ = V_series(0.1, 2, h=arb(1)/32, T=4.6)
vraw = tompf(vf[0],60)
vc = Vcont_series(0.1, 1)
vcont01 = tompf(vc.coeffs()[0],60)
print("raw V(0.1) =", mp.nstr(vraw,12))
print("cont V(0.1)=", mp.nstr(vcont01,12))
