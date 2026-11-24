from mpmath import mp, mpf, mpc, quad, sqrt, exp, log, pi, cos, sin
mp.dps = 150
p = 1 + mpc(0, mpf('0.12')); e = 1/p
g = lambda u: exp(-u)*sqrt((e*exp(-(e-1)*u))**2 + (e*(1-exp(-u))**(e-1))**2)
bp1 = [log(2), 2, 4, 6, 8, 12, 16, 24, 40, 80, mp.inf]
r1 = 8*quad(g, bp1, method='tanh-sinh')
r2 = 8*quad(g, bp1, method='tanh-sinh', degree=12)
bp2 = [log(2)+mpf('0.3'), 2.5, 5, 7.5, 10, 14, 20, 30, 60, 120, mp.inf]
r3 = 8*quad(g, bp2, method='tanh-sinh', degree=10)
print("r1:", mp.nstr(r1, 50))
print("r2:", mp.nstr(r2, 50))
print("r3:", mp.nstr(r3, 50))
print("spread:", mp.nstr(max(abs(r1-r2), abs(r1-r3), abs(r2-r3)), 3))
