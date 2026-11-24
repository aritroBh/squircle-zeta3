from mpmath import mp, mpf, pi, log, pslq, nstr
mp.dps = 55
import json
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2)
# n=13 k=2 with big maxcoeff
v = mpf(res["13,2"])
r = pslq([v,mpf(1),ln2,pi**2], tol=mpf('1e-40'), maxcoeff=10**12, maxsteps=10000)
print("13,2:", r)
if r:
    print("resid:", nstr(r[0]*v+r[1]+r[2]*ln2+r[3]*pi**2,3))
# n=1 k=2 with several bases
v1 = mpf(res["1,2"])
for basis,names in [([mpf(1),ln2,pi**2],["1","ln2","pi^2"]),
                    ([mpf(1),ln2,ln2**2,pi**2],["1","ln2","ln2^2","pi^2"]),
                    ([mpf(1),ln2,ln2**2,pi**2,pi**2*ln2],["1","ln2","ln2^2","pi^2","pi^2 ln2"])]:
    r = pslq([v1]+basis, tol=mpf('1e-38'), maxcoeff=10**10, maxsteps=10000)
    print("1,2", names, r)
