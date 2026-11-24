import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, nstr
mp.dps = 80
def h_lame(theta, p, lam):
    """support function of |x|^p+|y|^{lam p}=1, theta in [0,pi/2]."""
    if theta == 0: return mpf(1)
    if theta == pi/2: return mpf(1)
    c, s = cos(theta), sin(theta)
    def hof(a):
        return exp(-a/p)*c + (1-exp(-a))**(1/(lam*p))*s
    def F(a):
        return -exp(-a/p)*c + (1/lam)*(1-exp(-a))**(1/(lam*p)-1)*exp(-a)*s
    # initial guess: asymptotics for small/large a
    guesses = [mpf('1e-6'), mpf('1e-3'), mpf('0.1'), mpf(1), mpf(5), mpf(20), mpf(60)]
    best = None
    for g in guesses:
        try:
            r = mpmath.findroot(F, g, tol=mpf(10)**(-60), verify=False)
            if r >= 0:
                hv = hof(r)
                if best is None or hv > best[1]: best = (r, hv)
        except Exception:
            pass
    if best is None:
        return max(c, s)
    return max(best[1], max(c, s))
def P_lame_sup(p, lam):
    return 4*mpmath.quad(lambda th: h_lame(th, p, lam), [0, pi/8, pi/4, 3*pi/8, pi/2])
def P_lame_arc(p, lam):
    a = 2/p; b = 2/(lam*p)
    def f(t):
        dx = -a*cos(t)**(a-1)*sin(t)
        dy = b*sin(t)**(b-1)*cos(t)
        return sqrt(dx**2+dy**2)
    return 4*mpmath.quad(f, [0, pi/4, pi/2])
# validate at moderate p
for p,lam in [(mpf(5),mpf(2)), (mpf(3),mpf('0.5')), (mpf(20),mpf(2))]:
    Pa = P_lame_sup(p,lam); Pb = P_lame_arc(p,lam)
    print(f"p={p}, lam={lam}: support {nstr(Pa,45)}")
    print(f"            arclength {nstr(Pb,45)}  diff {nstr(Pa-Pb,4)}")
