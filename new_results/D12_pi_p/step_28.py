
import mpmath as mp
import pickle
from math import gcd
mp.mp.dps = 95
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi=mp.pi; L2=mp.log(2); z3=mp.zeta(3)

def identify(tag, t, b, nm):
    rel = mp.pslq([t]+b, maxcoeff=10**10, maxsteps=10000)
    if not rel:
        print(f"{tag}: no relation"); return
    g=0
    for x in rel: g=gcd(g,x)
    rel=[x//g for x in rel]
    val = -sum(rel[i+1]*b[i] for i in range(len(b)))/rel[0]
    terms = " + ".join(f"{rel[i+1]}*{nm[i]}" for i in range(len(b)) if rel[i+1])
    print(f"{tag}: {rel[0]}*({tag}) + [{terms}] = 0 ; residual={mp.nstr(abs(val-t),2)}")

identify('N2', N[2], [L2**3, pi**2], ['log2^3','pi^2'])
identify('N3', N[3], [L2**4, pi**2, pi**2*L2, z3], ['log2^4','pi^2','pi^2 log2','zeta3'])
identify('N4', N[4], [L2**5, pi**2, pi**2*L2, pi**2*L2**2, z3, z3*L2, pi**4],
         ['log2^5','pi^2','pi^2 log2','pi^2 log2^2','zeta3','zeta3 log2','pi^4'])

print("\nsqrt2-absence test for N3 (basis: Euclidean constants + ours):")
sq2 = mp.sqrt(2); las = mp.log(1+sq2); li2m = mp.polylog(2, sq2-1)
bE = [sq2, las, sq2*las, li2m, las**2, L2**4, pi**2, pi**2*L2, z3]
rel = mp.pslq([N[3]]+bE, maxcoeff=10**9, maxsteps=10000)
if rel:
    g=0
    for x in rel: g=gcd(g,x)
    rel=[x//g for x in rel]
    print("relation:", rel)
    print("-> uses sqrt2 constants? coefficients for [sq2,las,sq2*las,li2m,las^2]:", rel[1:6])
