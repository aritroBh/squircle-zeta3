code = r'''
from squircle import *
from methodC import *
from c1 import c1_direct, sigma_exact
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum
import json
setdps(65)
G = {}
G.update(json.load(open("taugrid34.json")))
G.update(json.load(open("taugrid56.json")))
def getZ(tau, s_): return mpf(G[f"{tau},{s_}"])
s_ = 4
for dt in [0.03, 0.06, 0.1, 0.15]:
    fd = (getZ(dt, s_) - getZ(-dt, s_))/(2*dt)
    print(f"FD dt={dt}: {nstr(fd, 22)}")
print("c1 lattice:", nstr(c1_direct(s_, K=800), 22))
# brute c1 lattice sum to cross-check c1_direct
K = 3000
P = fsum([power(k, -s_-1)*sigma_exact(k) for k in range(2, K+1)])
# tail via EM: sigma(k) ~ (k^2/2)ln k - k^2/4 - (k ln k)/2 + ... use integral approx for the rest
import mpmath
def sigk(k):
    n = k-1
    return (n*n/2+n/2+mpf(1)/12)*log(n) - n*n/4 + log(mpmath.glaisher)
tail = fsum([power(k, -s_-1)*sigk(k) for k in range(K+1, 30*K+1)])
print("brute c1 (K=3000 + asym tail to 90000):", nstr(-8*s_*(P+tail), 22))
'''
exec(code)
