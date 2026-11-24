from mpmath import mp, mpf, mpc, quad, sqrt, exp, log, pi
mp.dps = 110
p = 1 + mpc(0, mpf('0.12')); e = 1/p
f = lambda u: exp(-u)*sqrt((e*exp(-(e-1)*u))**2 + (e*(1-exp(-u))**(e-1))**2)
bps = [log(2), 2, 4, 6, 8, 12, 16, 24, 40, 80, 160, mp.inf]
for md in [8, 10, 12, 14]:
    try:
        r = 8*quad(f, bps, method='tanh-sinh', maxdegree=md)
        print("maxdegree", md, mp.nstr(r, 45))
    except Exception as ex:
        print("maxdegree", md, "failed:", ex)
