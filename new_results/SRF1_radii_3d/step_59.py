
import sys, time
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, arb_series, ctx
from vlib import ts_nodes, ser, tompf, V_series
flint.ctx.cap = 32
ctx.dps = 50; mp.mp.dps = 50

def expm1_ser(y, nlen):
    R = arb_series([1], prec=nlen)
    for k in range(nlen+2, 1, -1):
        R = 1 + y*R/k
    return y*R

def Vcont_series(eps0, nlen=2, h=arb(1)/32, T=4.2):
    X = ser(eps0, nlen)
    nodes = ts_nodes(h, T)
    cpre = [(c, c*c*(X*(-2*c.log())).exp(), (X*c.log()).exp()/c, wc) for (c,cm,wc,pc) in nodes]
    full = arb_series([0], prec=nlen)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(X*(-2*la)).exp()
        S0 = (1+A2).sqrt()
        B1 = -(1+2*X)*a*(X*(-2*l1a)).exp()/(1+a)**2*S0
        prefa = (X*(2*la)).exp()/a
        row = arb_series([0], prec=nlen)
        for (c, C2, prefc, wc) in cpre:
            u = a*c/(1+a)
            S = (1+A2+A2*C2).sqrt()
            T1 = expm1_ser((-1-2*X)*u.log1p(), nlen)*S
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

vf, _ = V_series(0.1, 2, h=arb(1)/32, T=4.6)
print("raw V(0.1) =", mp.nstr(tompf(vf[0],60),12))
print("cont V(0.1)=", mp.nstr(tompf(Vcont_series(0.1, 1).coeffs()[0],60),12))
print("(e+1)*V(e) as e -> -1:")
for e in [-0.75, -0.85, -0.9, -0.95]:
    v0 = tompf(Vcont_series(e, 1).coeffs()[0],60)
    print(f"e={e}: V={mp.nstr(v0,10)}  (e+1)V={mp.nstr((e+1)*v0,10)}")
