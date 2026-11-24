
import mpmath as mp, flint
mp.mp.dps = 100
alpha = mp.pi/3; c = mp.cos(alpha)
K = 300
flint.ctx.prec = 768; flint.ctx.cap = K+1
ca = flint.arb(mp.nstr(c, 105))
s = flint.acb_series([0,1], prec=K+1)

for xval in ['1.5', '30', '300', '1200']:
    x_m = mp.mpf(xval)
    ex_m = mp.exp(-x_m); L_m = mp.log1p(ex_m)
    Bt0_m = (ex_m-2*c)/(1+mp.sqrt(1-2*c*ex_m+ex_m**2))
    g0_m = Bt0_m/(mp.exp(x_m)+1)   # k=0 coefficient
    x = flint.acb(mp.nstr(x_m, 105))
    ex = (-x).exp(); L = (1+ex).log()
    eL = ((-L)*s).exp()
    Z = ex*(x*s).exp()
    SQ = (1 - 2*ca*Z + Z*Z).sqrt()
    Bt = (Z - 2*ca)/(1 + SQ)
    G = eL*Bt/(x.exp()+1)
    cf = G.coeffs()
    g0_f = cf[0]
    print("x=%-5s k0: flint=%s  err=%.2e  ball=%s" % (xval, g0_f.mid().real.str(25), abs(g0_f.mid().real-mp.mpf(mp.nstr(g0_m,40))), g0_f.rad().str(3)))
