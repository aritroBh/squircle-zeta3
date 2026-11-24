"""
H1, octahedron end.  Direct numerical surface area of |x|^p+|y|^p+|z|^p = 1, independent of the
Dirichlet/Gamma machinery.

Method N1 (graph of z over the piece where z is the largest coordinate):
  by the 48-fold symmetry (8 octants x 3! orderings) S(p) = 48 * Int_R sqrt(1 + (x/z)^{2(p-1)} + (y/z)^{2(p-1)}) dx dy,
  R = {0 <= y <= x, z >= x},  z = (1 - x^p - y^p)^{1/p}   (z_x = -(x/z)^{p-1}, z_y = -(y/z)^{p-1}).
  R = {0<=x<=3^{-1/p}, 0<=y<=x}  U  {3^{-1/p}<=x<=2^{-1/p}, 0<=y<=(1-2x^p)^{1/p}};
  inner variable y = (upper limit)*s, s in [0,1]; nested tanh-sinh.
Method N2 (projection of one octant onto the xy-plane, u = x^p, v = y^p, w = 1-u-v):
  S(p) = 8 eps^2 Int_0^1 (1-u) Int_0^1 (u v w)^{-eta} sqrt(u^{2eta}+v^{2eta}+w^{2eta}) dt du,
  v = (1-u) t, eps = 1/p, eta = 1 - eps.
Usage: python 03_direct_quadrature.py dps  (writes direct_S_dps{dps}.json)
"""
import sys, json
from multiprocessing import Pool
import mpmath as mp


def S_N1(p, dps):
    mp.mp.dps = dps
    p = mp.mpf(p)
    q = 2 * (p - 1)

    def Fw(x, y, w):
        # w = z^p = 1 - x^p - y^p, passed in a cancellation-free form
        lz = mp.log(w) / p
        return mp.sqrt(1 + mp.exp(q * (mp.log(x) - lz)) + mp.exp(q * (mp.log(y) - lz)))

    xs = mp.mpf(3) ** (-1 / p)
    xe = mp.mpf(2) ** (-1 / p)

    def g1(x):
        xp = x ** p
        return x * mp.quad(lambda s: Fw(x, x * s, 1 - xp * (1 + s ** p)), [0, 1])

    def g2(x):
        xp = x ** p
        Yp = max(1 - 2 * xp, mp.mpf(0))
        if Yp == 0:
            return mp.mpf(0)
        Y = Yp ** (1 / p)
        return Y * mp.quad(lambda s: Fw(x, Y * s, xp + Yp * (1 - s ** p)), [0, 1])

    I1 = mp.quad(g1, [0, xs])
    I2 = mp.quad(g2, [xs, xe])
    return 48 * (I1 + I2)


def S_N2(p, dps):
    mp.mp.dps = dps
    p = mp.mpf(p)
    eps = 1 / p
    eta = 1 - eps

    def inner(u):
        def h(t):
            v = (1 - u) * t
            w = (1 - u) * (1 - t)
            return (u * v * w) ** (-eta) * mp.sqrt(u ** (2 * eta) + v ** (2 * eta) + w ** (2 * eta))
        return (1 - u) * mp.quad(h, [0, 1])

    return 8 * eps ** 2 * mp.quad(inner, [0, 1])


def job(args):
    dstr, dps, meth = args
    mp.mp.dps = dps
    p = 1 + mp.mpf(dstr)
    val = S_N1(p, dps) if meth == 'N1' else S_N2(p, dps)
    return dstr, meth, mp.nstr(val, dps)


if __name__ == '__main__':
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    ds = ['0', '1', '-0.1', '0.1', '-0.03', '0.03', '-0.01', '0.01', '-0.003', '0.003', '-0.001', '0.001']
    jobs = [(d, dps, m) for d in ds for m in ('N1', 'N2') if not (d == '1' and m == 'N2')]  # N2 is slow at p=2
    with Pool(14) as pool:
        res = pool.map(job, jobs)
    out = {}
    for dstr, meth, val in res:
        out.setdefault(dstr, {})[meth] = val
    mp.mp.dps = dps
    print('exact checks: 4*sqrt(3) =', mp.nstr(4 * mp.sqrt(3), dps), ' 4*pi =', mp.nstr(4 * mp.pi, dps))
    for dstr in ds:
        a = mp.mpf(out[dstr]['N1'])
        b = mp.mpf(out[dstr]['N2']) if 'N2' in out[dstr] else mp.nan
        print(f'd={dstr:>7}  N1={out[dstr]["N1"]}  |N1-N2|={mp.nstr(abs(a - b), 3)}')
    json.dump(out, open(f'direct_S_dps{dps}.json', 'w'), indent=1)
    print('saved', f'direct_S_dps{dps}.json')
