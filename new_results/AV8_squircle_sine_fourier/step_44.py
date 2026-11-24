from mpmath import mp, mpf, pi, log, zeta, pslq, catalan, euler, nstr
mp.dps = 60
import json
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2); G = catalan
beta3 = mpf(0)  # Dirichlet beta(3) = sum (-1)^k/(2k+1)^3
from mpmath import lerchphi, inf, nsum
beta3 = nsum(lambda k: (-1)**k/(2*k+1)**3, [0, inf])
print("beta(3) =", nstr(beta3,20), " pi^3/32 =", nstr(pi**3/32,20))

ns=[1,3,5,7,9,11,13]
print("\n=== k=1 ===")
for n in ns:
    v = mpf(res[f"{n},1"])
    r = pslq([v, 1, ln2], tol=mpf('1e-45'), maxcoeff=10**8)
    print(n, nstr(v,25), r)
