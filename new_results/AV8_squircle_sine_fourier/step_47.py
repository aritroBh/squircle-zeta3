from mpmath import mp, mpf, pi, log, zeta, pslq, catalan, nstr
mp.dps = 55
import json
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2); G = catalan
ns=[1,5,7,9,11,13]
for tol in [mpf('1e-40'), mpf('1e-35')]:
    print(f"--- tol={nstr(tol,2)} ---")
    for n in ns:
        v = mpf(res[f"{n},2"])
        r = pslq([v,mpf(1),ln2,ln2**2,pi**2,pi**2*ln2], tol=tol, maxcoeff=10**9)
        if r:
            resid = r[0]*v + sum(a*b for a,b in zip(r[1:],[mpf(1),ln2,ln2**2,pi**2,pi**2*ln2]))
            print(n, r, "resid:", nstr(resid,3))
        else:
            print(n, "none")
