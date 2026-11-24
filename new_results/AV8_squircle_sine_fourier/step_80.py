from mpmath import mp, mpf, pi, log, zeta, pslq, catalan, nstr
mp.dps = 55
import json, itertools
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2); G = catalan
v3 = mpf(res["3,3"])
C = {"1":mpf(1),"ln2":ln2,"ln2^2":ln2**2,"ln2^3":ln2**3,"pi^2":pi**2,"pi^2 ln2":pi**2*ln2,
     "pi^2 ln2^2":pi**2*ln2**2,"zeta3":zeta(3),"zeta3 ln2":zeta(3)*ln2,"pi^3":pi**3}
best=[]
for sz in [4,5,6,7]:
    for combo in itertools.combinations(C.keys(), sz):
        basis=[C[c] for c in combo]
        r = pslq([v3]+basis, tol=mpf('1e-45'), maxcoeff=10**8, maxsteps=5000)
        if r:
            resid = r[0]*v3+sum(a*b for a,b in zip(r[1:],basis))
            if abs(resid)<mpf('1e-47'):
                best.append((abs(resid), combo, r))
for resid,combo,r in sorted(best)[:15]:
    print(combo, r, nstr(resid,3))
