from mpmath import mp, mpf, pi, log, zeta, pslq, catalan, nstr
mp.dps = 70
import json
res = json.load(open('scratchpad/bnder_hi.json'))
ln2 = log(2); G = catalan
v = mpf(res["3,3"])
C = {"1/pi":1/pi, "zeta3/pi":zeta(3)/pi, "1":mpf(1), "ln2":ln2, "ln2^2":ln2**2, "ln2^3":ln2**3,
     "pi":pi, "pi^2":pi**2, "pi^2 ln2":pi**2*ln2, "pi^3":pi**3, "zeta3":zeta(3), "zeta3 ln2":zeta(3)*ln2,
     "G":G, "pi G":pi*G}
import itertools
hits=[]
for sz in [7,8,9]:
    for combo in itertools.combinations(C.keys(), sz):
        basis=[C[c] for c in combo]
        r = pslq([v]+basis, tol=mpf('1e-60'), maxcoeff=10**7, maxsteps=5000)
        if r:
            resid = r[0]*v+sum(a*b for a,b in zip(r[1:],basis))
            if abs(resid) < mpf('1e-62'):
                hits.append((combo,r,nstr(resid,3)))
for h in hits[:10]: print(h)
print("done, hits:", len(hits))
