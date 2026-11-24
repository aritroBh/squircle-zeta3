from mpmath import mp, mpf, pi, log, pslq, nstr
mp.dps = 55
import json
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2)
for key in ["15,2","17,2"]:
    v = mpf(res[key])
    r = pslq([v,mpf(1),ln2,pi**2], tol=mpf('1e-40'), maxcoeff=10**12, maxsteps=10000)
    print(key, r)
    if r:
        print("   resid:", nstr(r[0]*v+r[1]+r[2]*ln2+r[3]*pi**2,3))
