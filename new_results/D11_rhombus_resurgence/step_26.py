from mpmath import mp, mpf, mpc, quad, sqrt, exp, log, pi
mp.dps = 100
p = 1 + mpc(0, mpf('0.12')); e = 1/p
# model: int e^{-e u} du from log2 to inf = e^{-e log2}/e
m1 = quad(lambda u: exp(-e*u), [log(2), 2, 4, 6, 8, 12, 16, 24, 40, 80, mp.inf], method='tanh-sinh')
print("model TS:", mp.nstr(m1, 40), " exact:", mp.nstr(exp(-e*log(2))/e, 40), " err:", mp.nstr(abs(m1-exp(-e*log(2))/e), 3))
# where does |w| get small? scan
g = lambda u: (e*exp(-(e-1)*u))**2 + (e*(1-exp(-u))**(e-1))**2
import numpy as np
us = np.linspace(float(log(2)), 30, 2000)
ws = [g(mpf(u)) for u in us]
aw = [abs(w) for w in ws]
i0 = int(np.argmin(aw))
print("min |w| =", mp.nstr(aw[i0], 5), "at u =", us[i0], " arg w =", mp.nstr(mp.arg(ws[i0])/pi, 5), "pi")
# try real/imag separate quadrature of full integrand with tight error tolerance
mp.quad = mp.quad  # noop
f = lambda u: exp(-u)*sqrt((e*exp(-(e-1)*u))**2 + (e*(1-exp(-u))**(e-1))**2)
fr = quad(lambda u: f(u).real, [log(2), 2, 4, 6, 8, 12, 16, 24, 40, 80, 160, mp.inf], method='tanh-sinh', error=True)
fi = quad(lambda u: f(u).imag, [log(2), 2, 4, 6, 8, 12, 16, 24, 40, 80, 160, mp.inf], method='tanh-sinh', error=True)
print("re:", mp.nstr(fr[0], 40), " est err:", mp.nstr(fr[1], 3))
print("im:", mp.nstr(fi[0], 40), " est err:", mp.nstr(fi[1], 3))
