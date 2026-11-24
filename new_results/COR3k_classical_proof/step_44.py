import mpmath as mp
mp.mp.dps = 120
def bracket(x):
    lam = mp.log(2); lx = mp.log(x); l1m = mp.log(1-x); l1p = mp.log(1+x)
    return (-mp.polylog(2, x) + mp.polylog(2, -x) + mp.polylog(2, (1+x)/2)
            - mp.polylog(2, (1-x)/(1+x)) - lx*l1m + lx*l1p + l1m*l1p
            - l1m*lam - l1p**2/2 + lam**2/2 + mp.pi**2/12)
for xv in [mp.mpf('0.13'), mp.mpf('0.377'), mp.mpf('0.61'), mp.mpf('0.95')]:
    print(xv, mp.nstr(bracket(xv), 4))
