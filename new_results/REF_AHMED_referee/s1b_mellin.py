# Re-run of the Mellin identity of thm:P2 with the cancellation-free integrand k_c = sin^2(th)/(s+c+r).
# (s1_P2.py used r-s-c directly; its large-s cancellation made that block meaningless.)
from mpmath import mp, mpf, quad, sqrt, cos, sin, pi, hyp2f1, inf, nstr, fabs
mp.dps = 50
out = []
def Lam(eps, th):
    z = sin(th/2)**2; return 2*z*hyp2f1(-eps, 1+eps, 2, z)
worst = 0
for th in [mpf('0.3'), pi/3, pi/2, 2*pi/3, 3*pi/4, mpf('2.9'), mpf('3.1')]:
    c = cos(th)
    for eps in [mpf('0.05'), mpf('0.37'), mpf('0.5'), mpf('0.83'), mpf('0.97')]:
        f = lambda s: sin(th)**2/(s + c + sqrt(1 + 2*c*s + s*s)) * s**(eps-1)
        g = lambda t: sin(th)**2/(1 + c*t + sqrt(t*t + 2*c*t + 1)) * t**(-eps)   # s=1/t on [1,inf)
        k0 = 1 - c; g0 = sin(th)**2/2   # endpoint values, integrated exactly
        a = (quad(lambda s: f(s) - k0*s**(eps-1), [0, 1]) + k0/eps
             + quad(lambda t: g(t) - g0*t**(-eps), [0, 1]) + g0/(1-eps))
        b = pi/sin(pi*eps)*Lam(eps, th)
        d = fabs(a-b)/fabs(b); worst = max(worst, d)
        out.append("th=%s eps=%s rel diff %s" % (nstr(th,6), nstr(eps,3), nstr(d,3)))
out.append("worst rel diff %s" % nstr(worst,3))
print("\n".join(out)); open('s1b_mellin.out','w').write("\n".join(out)+"\n")
