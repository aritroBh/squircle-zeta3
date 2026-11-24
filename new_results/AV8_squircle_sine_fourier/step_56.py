from mpmath import mp, mpf, pi, log, nsum, inf, pslq, nstr
mp.dps = 50
# using the conjectured first-order law for n>=3
def bp(n):
    return mpf((-1)**((n-1)//2))/((n**2-1)*(2*((n-1)//4)+1))
S = nsum(lambda j: bp(2*j+1)**2, [2, inf])  # n=2j+1>=5? careful j starts so that n>=3
# n=3 corresponds j=1
S = nsum(lambda j: bp(2*j+1)**2, [1, inf])
print("Sigma =", nstr(S,40))
r = pslq([S, mpf(1), pi**2], tol=mpf('1e-40'), maxcoeff=10**8)
print(r)
N2 = log(2)**2 - 3*log(2) + 5*pi**2/24
b1p = (log(2)-1)/2
b1pp = N2/2 - b1p**2 - S
print("b1'' predicted =", nstr(b1pp,40))
print("b1'' measured  =", nstr(mpf(__import__('json').load(open('scratchpad/bnder60.json'))["1,2"]),40))
