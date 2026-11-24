
import sys, time
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
import flint
from flint import arb, arb_series, ctx
from vlib import ts_nodes, ser, coeffs, tompf
flint.ctx.cap = 32
ctx.dps = 50; mp.mp.dps = 50

def Vcont_series(eps0, nlen=2, h=arb(1)/32, T=4.2):
    """V(eps) = Vsub + Pa + Pc near eps0 (Re eps0 > -1). arb_series version."""
    X = ser(eps0, nlen)
    nodes = ts_nodes(h, T)
    # precompute c-side
    cpre = []
    for (c, cm, wc, pc) in nodes:
        lc = c.log()
        C2 = c*c*(X*(-2*lc)).exp()
        prefc = (X*lc).exp()/c            # c^{e-1}
        A1c = -(1+2*X)*((X*(-(2+2*X.numerator if False else 0))) )  # placeholder
        cpre.append((c, C2, prefc, wc, pc, lc))
    full = arb_series([0], prec=nlen)
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(X*(-2*la)).exp()
        g0 = (X*(-2*l1a)).exp()/(1+a)*(1+A2).sqrt()     # g(a,0)
        B1 = -(1+2*X)*a*(X*((-2-2*X).mid() if False else 0))  # placeholder
        prefa = (X*(2*la)).exp()/a                       # a^{2e-1}
        row = arb_series([0], prec=nlen)
        for (c, C2, prefc, wc, pc, lc) in cpre:
            N = 1 + a + a*c
            g = (X*(-2*N.log())).exp()/N*(1+A2+A2*C2).sqrt()
            A1 = -(1+2*X)*((X*(-(2)* (1+c).log())).exp()*(1+c)**-2 - 1)  # -(1+2e)[(1+c)^{-2-2e}-1]
            B1 = -(1+2*X)*a*(X*((-2)*l1a)).exp()/(1+a)**2*(1+A2).sqrt()
            row += prefc*(g - g0 - a*A1 - c*B1)*wc
        full += prefa*row*wa
    # Pa, Pc (1D integrals)
    Pa = arb_series([0], prec=nlen); Pc = arb_series([0], prec=nlen)
    for (c, cm, wc, pc, lc) in cpre:
        Ja = (X*lc).exp()/c*((X*(-2*(1+c).log())).exp()/(1+c)**2 - 1)
        Pa += -(1+2*X)/(2*X+1)*Ja*wc
    for (a, am, wa, pa) in nodes:
        la = a.log(); l1a = (1+a).log()
        A2 = a*a*(X*(-2*la)).exp()
        Ic = (X*(2*la)).exp()*(X*((-2)*l1a)).exp()/(1+a)**2*(1+A2).sqrt()
        Pc += -(1+2*X)/(X+1)*Ic*wa
    return full + Pa + Pc

print("(e+1)*V(e) as e -> -1 (arb, h=1/32):")
for e in [-0.75, -0.85, -0.9, -0.95]:
    v = Vcont_series(e, 2)
    cs = v.coeffs(); v0 = tompf(cs[0], 60)
    print(f"e={e}: V={mp.nstr(v0,10)}  (e+1)V={mp.nstr((e+1)*v0,10)}")
