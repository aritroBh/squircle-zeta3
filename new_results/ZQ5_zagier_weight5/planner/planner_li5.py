# Planner Li5 evaluator for real x, |x| < 1, at any precision (independent of ZQ5's w5lib):
#   |x| <= 1/2      : direct series;
#   1/2 < x < 1     : Li5(e^t) = sum_{k != 4} zeta(5-k) t^k/k! + (H_4 - log(-t)) t^4/4!,  t = log x  (|t| < 2 pi);
#   -1 < x < -1/2   : Li5(x) = Li5(x^2)/16 - Li5(-x).
import mpmath as mp
_zcache = {}
def _zeta_int(n):
    if n not in _zcache or _zcache[n][0] < mp.mp.prec:
        _zcache[n] = (mp.mp.prec, mp.zeta(n) if n != 1 else None)
    return _zcache[n][1]
def li5_series(x, eps):
    tot = mp.mpf(0); p = x; n = 1
    while abs(p) > eps:
        tot += p / mp.mpf(n)**5; n += 1; p *= x
    return tot
def li5_log(x, eps):
    t = mp.log(x); tot = mp.mpf(0); tk = mp.mpf(1); fact = mp.mpf(1); k = 0
    H4 = mp.mpf(25) / 12
    while True:
        if k == 4:
            term = (H4 - mp.log(-t)) * tk / fact
        else:
            z = _zeta_int(5 - k)
            term = z * tk / fact
        tot += term
        # zeta(5-k) = 0 for odd k >= 7; for even k the terms decrease like (|t|/2pi)^k, so stop at a small even-k term
        if k > 10 and k % 2 == 0 and abs(term) < eps: break
        k += 1; tk *= t; fact *= k
    return tot
def li5(x, digits):
    eps = mp.mpf(10)**(-(digits + 10))
    if abs(x) <= 0.5: return li5_series(x, eps)
    if x > 0: return li5_log(x, eps)
    return li5(x * x, digits) / 16 - li5(-x, digits)
if __name__ == '__main__':
    mp.mp.dps = 330
    for x in ('0.6', '0.9', '0.97', '-0.7', '-0.95'):
        X = mp.mpf(x)
        a = li5(X, 300); b = li5_series(X, mp.mpf(10)**(-310)); c = mp.polylog(5, X)
        print(x, 'log-exp vs series', mp.nstr(abs(a - b), 3), ' vs mp.polylog', mp.nstr(abs(a - c), 3))
    for x in ('0.9999', '0.99', '-0.9999'):
        X = mp.mpf(x)
        print(x, 'log-exp vs mp.polylog', mp.nstr(abs(li5(X, 300) - mp.polylog(5, X)), 3))
    # non-circular duplication check: Li5(y^2) = 16 (Li5(y) + Li5(-y)) with y = 0.7 (series for -0.7 directly)
    y = mp.mpf('0.7'); eps = mp.mpf(10)**(-310)
    print('dup y=0.7:', mp.nstr(li5(y*y, 300) - 16*(li5(y, 300) + li5_series(-y, eps)), 3))
