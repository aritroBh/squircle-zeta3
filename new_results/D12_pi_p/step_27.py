
import mpmath as mp
import pickle
mp.mp.dps = 80
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
pi=mp.pi; L2=mp.log(2); z3=mp.zeta(3)

print("=== Radius analysis ===")
print("k : |N_k|^(1/k)      k/(e*pi)      |N_{k+2}/N_k|^{1/2} * pi/(k+1)")
for k in range(10,41,2):
    rt = abs(N[k])**(1/k)
    q = mp.sqrt(abs(N[k+2]/N[k]))*pi/(k+1)
    print(f"{k:2d}  {mp.nstr(rt,10):>12}  {mp.nstr(k/(mp.e*pi),10):>12}   {mp.nstr(q,10)}")

print("\n=== PSLQ confirmations of closed forms (residuals at 90 digits) ===")
mp.mp.dps = 95
forms = {
 'N2': (N[2], [L2**3, pi**2], ['log2^3','pi^2']),
 'N3': (N[3], [L2**4, pi**2, pi**2*L2, z3], ['log2^4','pi^2','pi^2 log2','zeta3']),
 'N4': (N[4], [L2**5, pi**2, pi**2*L2, pi**2*L2**2, z3, z3*L2, pi**4], ['log2^5','pi^2','pi^2 log2','pi^2 log2^2','zeta3','zeta3 log2','pi^4']),
}
for tag,(t,b,nm) in forms.items():
    rel = mp.pslq([t]+b, tol=mp.mpf('1e-85'), maxcoeff=10**9, maxsteps=10000)
    from math import gcd
    g=0
    for x in rel: g=gcd(g,x)
    rel=[x//g for x in rel]
    val = -sum(rel[i+1]*b[i] for i in range(len(b)))/rel[0]
    print(f"{tag}: coeffs {rel}  residual={mp.nstr(abs(val-t),2)}")

print("\n=== sqrt2-absence test: is N3 in the sqrt2/Euclidean algebra? ===")
sq2 = mp.sqrt(2); las = mp.log(1+sq2); li2m = mp.polylog(2, sq2-1)
bE = [mp.mpf(1), sq2, las, sq2*las, li2m, las**2, L2, pi**2, z3]
rel = mp.pslq([N[3]]+bE, tol=mp.mpf('1e-80'), maxcoeff=10**8, maxsteps=10000)
print("PSLQ with Euclidean constants:", rel, "-> relation exists?" , rel is not None)
if rel:
    val = -sum(rel[i+1]*bE[i] for i in range(len(bE)))/rel[0]
    print("   residual:", mp.nstr(abs(val-N[3]),2), " max|coeff|:", max(abs(x) for x in rel))
