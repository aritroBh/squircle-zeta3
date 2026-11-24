from mpmath import mp, mpf, mpc, quad, sqrt, exp, log
mp.dps = 100
p = 1 + mpc(0, mpf('0.09')); e = 1/p
# t = exp(-u), u in [log 2, inf): integrand smooth, e^{-u} decay
g = lambda u: exp(-u)*sqrt((e*exp(-(e-1)*u))**2 + (e*(1-exp(-u))**(e-1))**2)
r_true = 8*quad(g, [log(2), 4, 10, 30, mp.inf], method='tanh-sinh')
print("u-variable TS:", mp.nstr(r_true, 55))
f = lambda t: sqrt((e*t**(e-1))**2 + (e*(1-t)**(e-1))**2)
r_ts = 8*quad(f, [0, mpf(1)/16, mpf(1)/4, mpf(1)/2], method='tanh-sinh')
r_gl = 8*quad(f, [0, mpf(1)/16, mpf(1)/4, mpf(1)/2], method='gauss-legendre', degree=12)
print("t TS err:", mp.nstr(abs(r_ts-r_true), 3), "  t GL err:", mp.nstr(abs(r_gl-r_true), 3))
# also check u-variable at two precisions for self-consistency
mp.dps = 140
g2 = lambda u: exp(-u)*sqrt((e*exp(-(e-1)*u))**2 + (e*(1-exp(-u))**(e-1))**2)
r_true2 = 8*quad(g2, [log(2), 4, 10, 30, mp.inf], method='tanh-sinh')
mp.dps = 100
print("u-var self-consistency 100 vs 140:", mp.nstr(abs(mpc(r_true)-mpc(r_true2)), 3))
