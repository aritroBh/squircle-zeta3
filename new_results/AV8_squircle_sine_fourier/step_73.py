from mpmath import mp, mpf, pi, log, zeta, pslq, nstr
mp.dps = 55
import json
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2)
basis = [mpf(1), ln2, ln2**2, ln2**3, pi**2, pi**2*ln2, zeta(3)]
names = ["1","ln2","ln2^2","ln2^3","pi^2","pi^2 ln2","zeta3"]
ns=[1,3,5,7,9,11,13]
for n in ns:
    v = mpf(res[f"{n},3"])
    r = pslq([v]+basis, tol=mpf('1e-44'), maxcoeff=10**8, maxsteps=10000)
    if r is None:
        print(f"n={n}: none")
    else:
        resid = r[0]*v + sum(a*b for a,b in zip(r[1:],basis))
        print(f"n={n}: lead={r[0]}, resid={nstr(resid,4)}, rest={r[1:]}")
