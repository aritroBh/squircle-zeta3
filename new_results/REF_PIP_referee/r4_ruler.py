"""Referee script 4: prop:ruler, fixed-q convergence.

L_q(p)/4 = 2 - eps F_q(eps),  eps = 1/p,
F_q(eps) = 2 int_0^{1/2} [ b - a((1+z)^{1/q} - 1) ] dt,  a=t^{eps-1}, b=(1-t)^{eps-1}, z=(t/(1-t))^{q(1-eps)}.
Taylor coefficients of F_q at 0 by trapezoid rule on |eps| = r (r<1), then the partial sums of the
1/p series are compared with the geometric value L_q(p) (x-form of r3) at p = 1.25, 1.5, 3.
"""
import mpmath as mp
mp.mp.dps = 34

def F(q, e):
    q = mp.mpf(q)
    def f(t):
        a = mp.power(t, e - 1); b = mp.power(1 - t, e - 1)
        z = mp.power(t / (1 - t), q * (1 - e))
        return b - a * mp.expm1(mp.log1p(z) / q)
    return 2 * mp.quad(f, [0, mp.mpf(1) / 64, mp.mpf(1) / 8, mp.mpf(1) / 2])

def xform(p, q):
    p, q = mp.mpf(p), mp.mpf(q)
    xs = mp.power(2, -1 / p)
    f = lambda x: (1 + (x**p / (1 - x**p))**(q * (p - 1) / p))**(1 / q)
    pts = [mp.mpf(0)] + [xs * (1 - mp.mpf(k) / (p * p * q)) for k in (4096, 512, 64, 8, 1)] + [xs]
    return 8 * mp.quad(f, sorted(set(x for x in pts if x >= 0)))

N = 400; r = mp.mpf('0.85')
for q in [1.5, 3, 4]:
    vals = [F(q, r * mp.expjpi(2 * mp.mpf(j) / N)) for j in range(N)]
    coef = []
    for k in range(N // 2):
        c = sum(vals[j] * mp.expjpi(-2 * mp.mpf(j * k) / N) for j in range(N)) / N / r**k
        coef.append(mp.re(c))
    # check real eps value
    print('q=%s  F(0.3) direct vs series: %s' % (q, mp.nstr(F(q, mp.mpf('0.3')) - sum(c * mp.mpf('0.3')**k for k, c in enumerate(coef)), 3)))
    print('   |c_k|^(1/k) at k=20,60,120,199:', [mp.nstr(abs(coef[k])**(mp.mpf(1) / k), 5) for k in (20, 60, 120, 199)])
    for p in [3, 1.5, 1.25]:
        e = 1 / mp.mpf(p)
        direct = xform(p, q)
        for K in (10, 40, 199):
            ser = 4 * (2 - e * sum(coef[k] * e**k for k in range(K + 1)))
            print('   p=%-5s K=%-3d series - L_q(p) = %s' % (p, K, mp.nstr(ser - direct, 3)))
