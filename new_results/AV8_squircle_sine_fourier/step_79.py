from mpmath import mp, mpf, pi, log, zeta, pslq, catalan, nstr
mp.dps = 55
import json
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2); G = catalan
v3 = mpf(res["3,3"])
bases = [
    (["1","ln2","ln2^2","ln2^3","pi^2","pi^2 ln2","pi^3","zeta3"],
     [mpf(1),ln2,ln2**2,ln2**3,pi**2,pi**2*ln2,pi**3,zeta(3)]),
    (["1","ln2","ln2^2","ln2^3","pi^2","pi^2 ln2","zeta3","zeta3 ln2"],
     [mpf(1),ln2,ln2**2,ln2**3,pi**2,pi**2*ln2,zeta(3),zeta(3)*ln2]),
    (["1","ln2","ln2^2","ln2^3","pi^2","pi^2 ln2","pi^2 ln2^2","zeta3"],
     [mpf(1),ln2,ln2**2,ln2**3,pi**2,pi**2*ln2,pi**2*ln2**2,zeta(3)]),
    (["1","ln2","ln2^2","ln2^3","pi^2","pi^2 ln2","pi^3","zeta3","G","pi G"],
     [mpf(1),ln2,ln2**2,ln2**3,pi**2,pi**2*ln2,pi**3,zeta(3),G,pi*G]),
]
for names,basis in bases:
    r = pslq([v3]+basis, tol=mpf('1e-42'), maxcoeff=10**7, maxsteps=10000)
    if r:
        resid = r[0]*v3+sum(a*b for a,b in zip(r[1:],basis))
        print(len(basis), names, "\n  ", r, "resid:", nstr(resid,4))
    else:
        print(len(basis), names, "none")
