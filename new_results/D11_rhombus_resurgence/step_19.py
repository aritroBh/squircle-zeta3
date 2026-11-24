from mpmath import mp, mpf, mpc, quad, sqrt
mp.dps = 100
p = 1 + mpc(0, mpf('0.09')); e = 1/p
f = lambda t: sqrt((e*t**(e-1))**2 + (e*(1-t)**(e-1))**2)
for deg in [4, 6, 8, 10]:
    r = 8*quad(f, [0, mpf(1)/16, mpf(1)/4, mpf(1)/2], method='tanh-sinh', degree=deg)
    print("TS deg", deg, mp.nstr(r, 55))
for deg in [8, 12, 16]:
    r = 8*quad(f, [0, mpf(1)/16, mpf(1)/4, mpf(1)/2], method='gauss-legendre', degree=deg)
    print("GL deg", deg, mp.nstr(r, 55))
