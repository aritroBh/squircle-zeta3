from mpmath import mp, mpf, pi, log, zeta, pslq, nstr
mp.dps = 50
ln2 = log(2)
S13 = mpf('-0.0584063758387517933882325979061622052448509')
names = ["1","ln2","ln2^2","pi^2","pi^2 ln2","pi^4"]
basis = [mpf(1),ln2,ln2**2,pi**2,pi**2*ln2,pi**4]
r = pslq([S13]+basis, tol=mpf('1e-42'), maxcoeff=10**10, maxsteps=10000)
print(r)
if r:
    print("resid:", nstr(r[0]*S13+sum(a*b for a,b in zip(r[1:],basis)),4))
