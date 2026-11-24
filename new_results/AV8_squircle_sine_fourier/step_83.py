from mpmath import mp, mpf, pi, log, zeta, pslq, catalan, nstr
mp.dps = 70
import json
res = json.load(open('scratchpad/bnder_hi.json'))
ln2 = log(2); G = catalan
names = ["1","ln2","ln2^2","ln2^3","pi^2","pi^2 ln2","pi^2 ln2^2","zeta3","zeta3 ln2","pi^3","G","pi G"]
basis = [mpf(1),ln2,ln2**2,ln2**3,pi**2,pi**2*ln2,pi**2*ln2**2,zeta(3),zeta(3)*ln2,pi**3,G,pi*G]
for key in ["1,3","3,3","5,3"]:
    v = mpf(res[key])
    r = pslq([v]+basis, tol=mpf('1e-58'), maxcoeff=10**7, maxsteps=10000)
    print(key, "->", r if r is None else (r, nstr(r[0]*v+sum(a*b for a,b in zip(r[1:],basis)),4)))
