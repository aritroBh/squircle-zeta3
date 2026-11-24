from mpmath import mp, mpf, pi, log, zeta, pslq, nstr
mp.dps = 50
ln2 = log(2)
def bp(n):
    n=int(n)
    if n==1: return (ln2-1)/2
    return mpf((-1)**((n-1)//2))/((n**2-1)*(2*((n-1)//4)+1))
bpp = {1: mpf(3)/8 - ln2 + ln2**2/4 + 5*pi**2/96 - pi**4/768,
       3: (12 - 16*ln2 + 3*pi**2)/128,
       5: (508 + 48*ln2 - 45*pi**2)/3456,
       7: (-572 - 96*ln2 + 63*pi**2)/13824,
       9: (2532 + 160*ln2 - 225*pi**2)/64000,
       11: (-25592 - 2160*ln2 + 2475*pi**2)/1296000,
       13: (1598344 + 75600*ln2 - 143325*pi**2)/88905600}
S13 = sum(bp(n)*bpp[n] for n in [1,3,5,7,9,11,13])
print("Sigma_13 =", nstr(S13,42))
names = ["1","ln2","ln2^2","pi^2","pi^2 ln2","pi^4","zeta3"]
basis = [mpf(1),ln2,ln2**2,pi**2,pi**2*ln2,pi**4,zeta(3)]
r = pslq([S13]+basis, tol=mpf('1e-40'), maxcoeff=10**9, maxsteps=10000)
print(r)
if r:
    print("resid:", nstr(r[0]*S13+sum(a*b for a,b in zip(r[1:],basis)),4))
