
import sys, time
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, arb_series, ctx
from vlib import ts_nodes, ser, tompf
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
        sq = (1+A2).sqrt()
        g0 = (X*(-2*l1a)).exp()/(1+a)*sq
        B1 = -(1+2*X)*a*(X*(-2*l1a)).exp()/(1+a)**2*sq
        prefa = (X*(2*la)).exp()/a
        row = arb_series([0], prec=nlen)
        for (c, C2, prefc, wc) in cpre:
            N = 1 + a + a*c
            g = (X*(-2*N.log())).exp()/N*(1+A2+A2*C2).sqrt()
            A1 = -(1+2*X)*c
            row += prefc*(g - g0 - a*A1 - c*B1)*wc
        full += prefa*row*wa
    Pa = -(1+2*X)/((2*X+1)*(X+1))
    Pc = arb_series([0], prec=nlen)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(X*(-2*la)).exp()
        Ic = (X*(2*la)).exp()*(X*(-2*l1a)).exp()/(1+a)**2*(1+A2).sqrt()
        Pc += -(1+2*X)/(X+1)*Ic*wa
    return full + Pa + Pc

print("(e+1)*V(e) as e -> -1 (corrected):")
for e in [-0.75, -0.85, -0.9, -0.95]:
    v = Vcont_series(e, 2)
    v0 = tompf(v.coeffs()[0], 60)
    print(f"e={e}: V={mp.nstr(v0,10)}  (e+1)V={mp.nstr((e+1)*v0,10)}")
# sanity at eps=0.1 vs raw vlib V_series
from vlib import V_series
vf, _ = V_series(0, 2, h=arb(1)/32, T=4.6)
V_raw_at = tompf(vf.coeffs()[0],60) + 0.1*tompf(vf.coeffs()[1],60)
vc = Vcont_series(0.1, 1)
print("sanity: raw V(0.1)~", mp.nstr(V_raw_at,10), "  Vcont(0.1)=", mp.nstr(tompf(vc.coeffs()[0],60),10))
