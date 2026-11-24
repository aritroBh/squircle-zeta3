from mpmath import mp, mpf, pi, log, zeta, pslq, nstr
mp.dps = 55
import json, itertools
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2)
consts = {"1":mpf(1),"ln2":ln2,"ln2^2":ln2**2,"pi^2":pi**2,"pi^2 ln2":pi**2*ln2}
for n in [1,5,7,9,11,13]:
    v = mpf(res[f"{n},2"])
    found=False
    for m in [2,3,4]:
        for combo in itertools.combinations(consts.keys(),m):
            basis=[consts[c] for c in combo]
            r = pslq([v]+basis, tol=mpf('1e-40'), maxcoeff=10**9)
            if r:
                resid = r[0]*v + sum(a*b for a,b in zip(r[1:],basis))
                if abs(resid) < mpf('1e-42'):
                    print(n, combo, r, "resid:", nstr(resid,3)); found=True
    if not found: print(n, "no small relation")
