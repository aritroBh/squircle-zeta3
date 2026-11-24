import time
import mpmath as mp
mp.mp.dps = 60

def hpl_coeffs(w, N):
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
            S += c[n - 1]                # includes c_0 (=1 only when rest=())
            out[n] = S / n
    else:
        A = mp.mpf(0)
        for n in range(1, N + 1):
            A = c[n - 1] - A             # A_{n-1} = sum_k (-1)^{n-1-k} c_k
            out[n] = A / n
    return out

def H_series(w, x, N=None):
    import math
    if N is None:
        N = int((mp.mp.dps + 15) / abs(math.log10(abs(x)))) + 10
    c = hpl_coeffs(w, N)
    tot = mp.mpf(0); xn = x
    for n in range(1, N + 1):
        tot += c[n] * xn; xn *= x
    return tot

x = mp.mpf('0.77')
print("H1 =", mp.nstr(H_series((1,), x) + mp.log(1-x), 4))
print("Hm1 =", mp.nstr(H_series((-1,), x) - mp.log(1+x), 4))
print("H01-Li2 =", mp.nstr(H_series((0,1), x) - mp.polylog(2, x), 4))
print("H0001-Li4 =", mp.nstr(H_series((0,0,0,1), x) - mp.polylog(4, x), 4))
print("H000m1+Li4(-x) =", mp.nstr(H_series((0,0,0,-1), x) + mp.polylog(4, -x), 4))
print("H0011 =", mp.nstr(H_series((0,0,1,1), x), 20))
print("H0101 =", mp.nstr(H_series((0,1,0,1), x), 20))
t0=time.time(); H_series((0,1,0,1), mp.mpf('0.966')); print(f"time at 0.966, dps60: {time.time()-t0:.2f}s")
