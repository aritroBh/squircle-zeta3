"""Referee script 3: independent geometric evaluation of pi_p and L_q(p).

x-parametrisation of the quarter arc up to the diagonal, y = (1-x^p)^{1/p}, |y'|^p = (x^p/(1-x^p))^{p-1}:
  pi_p   = 4 int_0^{2^{-1/p}} (1 + (x^p/(1-x^p))^{p-1})^{1/p} dx            (l^p length, half perimeter)
  L_q(p) = 8 int_0^{2^{-1/p}} (1 + (x^p/(1-x^p))^{q(p-1)/p})^{1/q} dx      (l^q length, full perimeter)
  L_inf(p) = 8 int_0^{2^{-1/p}} max(1,|y'|) dx = 8*2^{-1/p}  since |y'|<=1 there.
Compared with the paper's t-integrals and with the s-form split used in r1.
"""
import mpmath as mp
mp.mp.dps = 80

def xform(p, q):
    p, q = mp.mpf(p), mp.mpf(q)
    xs = mp.power(2, -1 / p)
    def f(x):
        r = x**p / (1 - x**p)
        return (1 + r**(q * (p - 1) / p))**(1 / q)
    # boundary layer near xs of width ~ 1/(p q)
    pts = [mp.mpf(0)] + [xs * (1 - mp.mpf(k) / (p * p * max(q, 1))) for k in (4096, 512, 64, 8, 1)] + [xs]
    pts = sorted(set([x for x in pts if x >= 0]))
    return 8 * mp.quad(f, pts)

def pi_x(p): return xform(p, p) / 2

def tform(p, q):  # paper's L_q, symmetric half, substitution t = y^p to remove t^{1/p-1}
    p, q = mp.mpf(p), mp.mpf(q)
    def f(y):
        t = y**p
        return (1 + (t / (1 - t))**(q * (1 - 1 / p)))**(1 / q)
    # (4/p) int_0^1 (...) dt = (8/p) int_0^{1/2} t^{1/p-1}(1+(t/(1-t))^{q(1-1/p)})^{1/q} dt, t=y^p: = 8 int_0^{2^{-1/p}} f(y) dy
    ys = mp.power(2, -1 / p)
    pts = [mp.mpf(0)] + [ys * (1 - mp.mpf(k) / (p * p * max(q, 1))) for k in (4096, 512, 64, 8, 1)] + [ys]
    pts = sorted(set([x for x in pts if x >= 0]))
    return 8 * mp.quad(f, pts)

def R_sform(p):
    p = mp.mpf(p); h = 1 / p; c = h / (1 - h)
    f = lambda s: mp.exp(-s * h * c) * (1 + mp.exp(-s * c)) ** (-1 - h) * mp.expm1(h * mp.log1p(mp.exp(-s)))
    return 4 * h * c * mp.quad(f, [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512])

print('pi_2 (x-form) - pi =', mp.nstr(pi_x(2) - mp.pi, 3))
for p in [3, 5, 10, 100, 1000, 10**4]:
    a = pi_x(p); b = 4 * mp.power(2, -1 / mp.mpf(p)) + R_sform(p)
    print('p=%-6s pi_p(x-form)=%s  split-form diff=%s' % (p, mp.nstr(a, 50), mp.nstr(a - b, 3)))
for p in [3, 5, 10]:
    p = mp.mpf(p); q = p / (p - 1)
    print('duality p=%s: pi_p - pi_q = %s' % (mp.nstr(p, 3), mp.nstr(pi_x(p) - pi_x(q), 3)))
print('min near 2:', [mp.nstr(pi_x(p), 20) for p in (1.9, 2, 2.1)])
# rulers
for p in [1, 1.5, 2, 3, 7, 50]:
    l1 = xform(p, 1) if p > 1 else mp.mpf(8)
    li = 8 * mp.power(2, -1 / mp.mpf(p))
    lbig = xform(p, 10**6) if p > 1 else None
    print('p=%-4s L_1 - 8 = %s ; L_q(q=1e6) - 8*2^{-1/p} = %s' % (p, mp.nstr(l1 - 8, 3) if p > 1 else 'exact', mp.nstr(lbig - li, 3) if lbig is not None else '-'))
def traw(p, q):  # paper's t-integral, raw (accurate only for small p)
    p, q = mp.mpf(p), mp.mpf(q)
    f = lambda t: (t**(q / p - q) + (1 - t)**(q / p - q))**(1 / q)
    return 4 / p * mp.quad(f, [0, mp.mpf(1)/4, mp.mpf(1)/2, mp.mpf(3)/4, 1])
for (p, q) in [(3, 2), (3, 3), (3, 1.5), (1.5, 4), (2, 3)]:
    print('L_%s(%s): x-form %s  raw t-integral diff %s' % (q, p, mp.nstr(xform(p, q), 40), mp.nstr(xform(p, q) - traw(p, q), 3)))
