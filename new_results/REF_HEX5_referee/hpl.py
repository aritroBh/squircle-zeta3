"""Referee's own harmonic polylogarithms H_w(x), letters 0 and 1, f_0 = 1/t, f_1 = 1/(1-t),
H_{a w}(x) = int_0^x f_a(t) H_w(t) dt, H_() = 1.  Only words ending in 1 (no trailing-zero regularisation).
Method: power series at 0 (|x| <= 1/2), then analytic continuation along the real axis by Taylor
re-expansion of the linear ODE system at successive centres (each step 0.4 of the distance to the
nearest singular point 0 or 1).  Valid for every real x < 1.
"""
import mpmath as mp
import itertools

WORDS = [w for n in range(1, 6) for w in itertools.product((0, 1), repeat=n) if w[-1] == 1]
WORDS.sort(key=len)


def _series0(nterms):
    """Taylor coefficients at 0 of all H_w (words ending in 1)."""
    h = {(): [mp.mpf(1)] + [mp.mpf(0)] * nterms}
    for w in WORDS:
        a, rest = w[0], w[1:]
        hr = h[rest]
        c = [mp.mpf(0)] * (nterms + 1)
        if a == 0:
            for n in range(1, nterms + 1):
                c[n] = hr[n] / n
        else:
            acc = mp.mpf(0)
            for n in range(0, nterms):
                acc += hr[n]
                c[n + 1] = acc / (n + 1)
        h[w] = c
    return h


def _evalser(c, y):
    s = mp.mpf(0)
    for a in reversed(c):
        s = s * y + a
    return s


def _step(vals, c0, y, nterms):
    """Given H_w(c0) for all words, return H_w(c0+y) via Taylor expansion at c0."""
    f0 = [mp.mpf((-1) ** k) / c0 ** (k + 1) for k in range(nterms + 1)]
    f1 = [1 / (1 - c0) ** (k + 1) for k in range(nterms + 1)]
    h = {(): [mp.mpf(1)] + [mp.mpf(0)] * nterms}
    out = {}
    for w in WORDS:
        a, rest = w[0], w[1:]
        hr = h[rest]
        f = f0 if a == 0 else f1
        c = [vals[w]] + [mp.mpf(0)] * nterms
        for n in range(0, nterms):
            c[n + 1] = mp.fsum(f[k] * hr[n - k] for k in range(n + 1)) / (n + 1)
        h[w] = c
        out[w] = _evalser(c, y)
    return out


def hpl_all(x):
    x = mp.mpf(x)
    assert x < 1
    nterms = int(3.5 * mp.mp.dps) + 20
    x0 = x if abs(x) <= mp.mpf(1) / 2 else (mp.mpf(1) / 2 if x > 0 else -mp.mpf(1) / 2)
    h = _series0(nterms)
    vals = {w: _evalser(h[w], x0) for w in WORDS}
    c0 = x0
    nt2 = int(2.6 * mp.mp.dps) + 20
    while c0 != x:
        rad = min(abs(c0), abs(1 - c0))
        stepmax = mp.mpf(2) / 5 * rad
        y = x - c0
        if abs(y) > stepmax:
            y = stepmax if y > 0 else -stepmax
        vals = _step(vals, c0, y, nt2)
        c0 = c0 + y
        if abs(x - c0) < mp.mpf(10) ** (-mp.mp.dps - 5):
            c0 = x
    return vals


ATOMS = {'H1': (1,), 'Li2': (0, 1), 'Li3': (0, 0, 1), 'H011': (0, 1, 1), 'Li4': (0, 0, 0, 1),
         'H0011': (0, 0, 1, 1), 'H0111': (0, 1, 1, 1), 'Li5': (0, 0, 0, 0, 1), 'H00011': (0, 0, 0, 1, 1),
         'H00101': (0, 0, 1, 0, 1), 'H00111': (0, 0, 1, 1, 1), 'H01011': (0, 1, 0, 1, 1), 'H01111': (0, 1, 1, 1, 1)}


def atoms(x):
    v = hpl_all(x)
    return {k: v[w] for k, w in ATOMS.items()}


if __name__ == '__main__':
    mp.mp.dps = 50
    for x in ['-0.3', '0.45', '-2', '-1', '0.8', '-7.5']:
        a = atoms(mp.mpf(x))
        xv = mp.mpf(x)
        errs = [abs(a['H1'] + mp.log(1 - xv))]
        for n, k in ((2, 'Li2'), (3, 'Li3'), (4, 'Li4'), (5, 'Li5')):
            errs.append(abs(a[k] - mp.polylog(n, xv)))
        # S_{1,2}(x) = H011 ; known: S12(x) = -Li3(1-x)+... use numeric integral instead
        s12 = mp.quad(lambda t: mp.log(1 - t) ** 2 / t, [0, xv]) / 2
        errs.append(abs(a['H011'] - s12))
        # H_{0011}: int_0^x H011(t)/t ... check H00101 against nested sum at |x|<1
        print(x, [mp.nstr(e, 3) for e in errs])
    # H00101 = Li_{3,2}(x,1) = sum_{n>m} x^n/(n^3 m^2): check at x=-0.3 by direct sum
    xv = mp.mpf('-0.3')
    s = mp.mpf(0); inner = mp.mpf(0)
    for n in range(1, 400):
        s += xv ** n * inner / mp.mpf(n) ** 3
        inner += mp.mpf(1) / n ** 2
    print('H00101(-0.3) vs nested sum', mp.nstr(atoms(xv)['H00101'] - s, 3))
    # x=-1 values vs nsum
    v = atoms(-1)
    print('Li5(-1)+15/16 z5', mp.nstr(v['Li5'] + mp.mpf(15) / 16 * mp.zeta(5), 3))
    print('H011(-1)-z3/8', mp.nstr(v['H011'] - mp.zeta(3) / 8, 3))
