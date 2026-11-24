"""hpl_series.py -- fast high-precision HPL evaluation via exact recursive series,
and multiple polylogs Li_{3,1}(x,y)=sum_{n>m} x^n y^m/(n^3 m), Li_{2,2} likewise.
Cross-checked against the quadrature of G1_zeta3_any_angle/hpl_eval.py (step15)."""
import math
import mpmath as mp

def hpl_coeffs(w, N):
    """coefficients c[0..N] of H_w(x) = sum_{n>=1} c_n x^n (c_0 = 0 for nonempty w)."""
    if len(w) == 0:
        c = [mp.mpf(0)] * (N + 1); c[0] = mp.mpf(1); return c
    a, rest = w[0], w[1:]
    c = hpl_coeffs(rest, N)
    out = [mp.mpf(0)] * (N + 1)
    if a == 0:
        for n in range(1, N + 1):
            out[n] = c[n] / n
    elif a == 1:
        S = mp.mpf(0)
        for n in range(1, N + 1):
            S += c[n - 1]
            out[n] = S / n
    else:
        A = mp.mpf(0)
        for n in range(1, N + 1):
            A = c[n - 1] - A
            out[n] = A / n
    return out

def H_series(w, x, N=None):
    w = tuple(w)
    if w == (0,):
        return mp.log(x)
    if N is None:
        N = int((mp.mp.dps + 15) / abs(math.log10(abs(x)))) + 10
    c = hpl_coeffs(w, N)
    tot = mp.mpf(0); xn = x
    for n in range(1, N + 1):
        tot += c[n] * xn; xn *= x
    return tot

def _mpl_series(x, y, a, b):
    tot = mp.mpf(0); inner = mp.mpf(0); xn = x; n = 1
    while True:
        if n >= 2:
            m = n - 1
            inner += y**m / mp.mpf(m)**b
        term = xn * inner / mp.mpf(n)**a
        tot += term
        if n > 40 and abs(term) * (1 + x / (1 - x) * 2) < mp.mpf(10) ** (-mp.mp.dps - 5):
            break
        n += 1; xn *= x
    return tot

def Li31(x, y): return _mpl_series(x, y, 3, 1)
def Li22(x, y): return _mpl_series(x, y, 2, 2)
