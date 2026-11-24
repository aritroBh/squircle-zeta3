from mpmath import mp, mpf, pi, log, zeta, pslq, nstr
mp.dps = 70
import json
res = json.load(open('scratchpad/bnder_hi.json'))
ln2 = log(2)
v = mpf(res["1,3"])
names = ["1","ln2","ln2^2","ln2^3","pi^2","pi^2 ln2","zeta3","pi^4"]
basis = [mpf(1),ln2,ln2**2,ln2**3,pi**2,pi**2*ln2,zeta(3),pi**4]
for tol in [mpf('1e-60'), mpf('1e-55')]:
    r = pslq([v]+basis, tol=tol, maxcoeff=10**7, maxsteps=10000)
    if r:
        resid = r[0]*v+sum(a*b for a,b in zip(r[1:],basis))
        print(tol, r, "resid:", nstr(resid,4))
    else:
        print(tol, "none")
