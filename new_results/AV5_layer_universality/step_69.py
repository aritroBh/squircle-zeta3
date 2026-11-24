import mpmath
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, nstr
mp.dps = 80
def h_brute(theta, p, lam, N=2000):
    c, s = cos(theta), sin(theta)
    # scan log-grid in a, then refine by golden section
    best = (-1, max(c,s))
    prev_a, prev_h = None, None
    import math
    xs = [exp(mpf(-20 + i*40/N)) for i in range(N+1)]
    hs = []
    for a in xs:
        hv = exp(-a/p)*c + (1-exp(-a))**(1/(lam*p))*s
        hs.append(hv)
    i0 = max(range(N+1), key=lambda i: hs[i])
    lo, hi = xs[max(0,i0-1)], xs[min(N,i0+1)]
    # golden section maximize
    gr = (sqrt(5)-1)/2
    x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
    f1 = exp(-x1/p)*c + (1-exp(-x1))**(1/(lam*p))*s
    f2 = exp(-x2/p)*c + (1-exp(-x2))**(1/(lam*p))*s
    for _ in range(300):
        if f1 < f2:
            lo, x1, f1 = x1, x2, f2
            x2 = lo+gr*(hi-lo)
            f2 = exp(-x2/p)*c + (1-exp(-x2))**(1/(lam*p))*s
        else:
            hi, x2, f2 = x2, x1, f1
            x1 = hi-gr*(hi-lo)
            f1 = exp(-x1/p)*c + (1-exp(-x1))**(1/(lam*p))*s
    return max(f1,f2)
def h_lame(theta, p, lam):
    c, s = cos(theta), sin(theta)
    def hof(a):
        return exp(-a/p)*c + (1-exp(-a))**(1/(lam*p))*s
    def F(a):
        return -exp(-a/p)*c + (1/lam)*(1-exp(-a))**(1/(lam*p)-1)*exp(-a)*s
    guesses = [mpf('1e-6'), mpf('1e-3'), mpf('0.1'), mpf(1), mpf(5), mpf(20), mpf(60)]
    best = max(c,s)
    for g in guesses:
        try:
            r = mpmath.findroot(F, g, tol=mpf(10)**(-60), verify=False)
            if r >= 0: best = max(best, hof(r))
        except Exception: pass
    return best
p, lam = mpf(5), mpf(2)
for th in [mpf('0.1'), mpf('0.4'), mpf('0.7853981633974483'), mpf('1.2'), mpf('1.5')]:
    print(nstr(th,8), nstr(h_lame(th,p,lam),40), nstr(h_brute(th,p,lam),40))
