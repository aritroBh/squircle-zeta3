from mpmath import mp, mpf, pi, log, zeta, pslq, nstr
mp.dps = 55
import json, itertools
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2)
C = {"1":mpf(1),"ln2":ln2,"ln2^2":ln2**2,"ln2^3":ln2**3,"pi^2":pi**2,"pi^2 ln2":pi**2*ln2,"pi^2 ln2^2":pi**2*ln2**2,"zeta3":zeta(3)}
# try increasingly large subsets, keep relations with residual < 1e-46 and small lead coeff
for n in [3]:
    v = mpf(res[f"{n},3"])
    print(f"=== n={n} value={nstr(v,30)} ===")
    for sz in [2,3,4,5,6,7]:
        hits=0
        for combo in itertools.combinations(C.keys(), sz):
            basis=[C[c] for c in combo]
            r = pslq([v]+basis, tol=mpf('1e-46'), maxcoeff=10**7)
            if r:
                resid = r[0]*v+sum(a*b for a,b in zip(r[1:],basis))
                if abs(resid)<mpf('1e-46') and max(abs(x) for x in r)<10**7:
                    print(sz, combo, r, nstr(resid,3)); hits+=1
                    if hits>3: break
        if hits: break
