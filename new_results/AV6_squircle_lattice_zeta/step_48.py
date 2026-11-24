code = r'''
from squircle import *
from mpmath import mp, mpf, nstr, log, power, fsum
setdps(50)
def inner(k, p):
    kf = mpf(k)
    return fsum([log(1+power(j/kf, p)) for j in range(1, k)])
for k in [100, 1000]:
    for p in [5, 10, 20]:
        v = inner(k,p)
        # integral approx
        I0q = fsum([(-1)**(n+1)/(n*(n*p+1)) for n in range(1, 200)])
        print(f"k={k} p={p}: inner={nstr(v,10)}  k*I0q={nstr(k*I0q,10)}  ratio={nstr(v/(k*I0q),6)}")
'''
exec(code)
