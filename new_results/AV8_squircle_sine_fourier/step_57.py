from mpmath import mp, mpf, pi, log, nsum, inf, pslq, nstr
mp.dps = 50
import json
def bp(n):
    n = int(n)
    return mpf((-1)**((n-1)//2))/((n**2-1)*(2*((n-1)//4)+1))
S = nsum(lambda j: bp(2*j+1)**2, [1, inf])
print("Sigma =", nstr(S,40))
r = pslq([S, mpf(1), pi**2], tol=mpf('1e-40'), maxcoeff=10**10)
print("pslq(S,1,pi^2):", r)
N2 = log(2)**2 - 3*log(2) + 5*pi**2/24
b1p = (log(2)-1)/2
b1pp = N2/2 - b1p**2 - S
res = json.load(open('scratchpad/bnder60.json'))
print("b1'' predicted =", nstr(b1pp,40))
print("b1'' measured  =", nstr(mpf(res["1,2"]),40))
print("diff:", nstr(b1pp - mpf(res["1,2"]),5))
