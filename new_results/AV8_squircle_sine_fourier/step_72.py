from mpmath import mp, mpf, pi, log, zeta, pslq, catalan, nstr
mp.dps = 55
import json
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2); G = catalan
# basis for k=3 (weight<=3 MZV-ish + level-4 controls)
basis = [mpf(1), ln2, ln2**2, ln2**3, pi**2, pi**2*ln2, pi**2*ln2**2, zeta(3), zeta(3)*ln2, pi**3, G, pi*G]
names = ["1","ln2","ln2^2","ln2^3","pi^2","pi^2 ln2","pi^2 ln2^2","zeta3","zeta3 ln2","pi^3","G","pi G"]
ns=[1,3,5,7,9,11,13]
for n in ns:
    v = mpf(res[f"{n},3"])
    r = pslq([v]+basis, tol=mpf('1e-40'), maxcoeff=10**9, maxsteps=10000)
    if r is None:
        print(f"n={n}: NO RELATION at 1e-40")
    else:
        resid = r[0]*v + sum(a*b for a,b in zip(r[1:],basis))
        print(f"n={n}: coeff {r[0]}, resid={nstr(resid,4)}")
        for a,nm in zip(r[1:],names):
            if a: print("     %+d*%s"%(a,nm))
