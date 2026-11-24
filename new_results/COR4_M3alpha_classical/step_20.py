import time, sys
import mpmath as mp
mp.mp.dps = 50

# Exact recursive series coefficients for H_w(x), letters {0,1,-1}, words not ending... all fine for |x|<1.
# H_w(x) = sum c_n x^n (n>=1).  Recursion on leftmost letter:
# H_{0,w}: c'_n = c_n / n
# H_{1,w}: c'_n = S_{n-1}/n where S_m = sum_{k<=m} c_k
# H_{-1,w}: c'_n = A_{n-1}/n where A_m = sum_{k<=m} (-1)^{k+1} c_k
def hpl_coeffs(w, N):
    # returns list c[0..N] of mpf coefficients of H_w
    if len(w) == 0:
        c = [mp.mpf(0)] * (N + 1)
        c[0] = mp.mpf(1)  # H_empty = 1
        return c
    a, rest = w[0], w[1:]
    c = hpl_coeffs(rest, N)
    out = [mp.mpf(0)] * (N + 1)
    if a == 0:
        for n in range(1, N + 1):
            out[n] = c[n] / n
    elif a == 1:
        S = mp.mpf(0)
        for n in range(1, N + 1):
            S += c[n - 1] if n - 1 >= 1 else 0
            # S = sum_{k=1}^{n-1} c_k
            out[n] = S / n
    else:
        A = mp.mpf(0)
        for n in range(1, N + 1):
            k = n - 1
            if k >= 1:
                A += ((-1) ** (k + 1)) * c[k]
            out[n] = A / n
    return out

def H_series(w, x, N=None):
    if N is None:
        # choose N so x^N < 1e-(dps+10)
        import math
        N = int((mp.mp.dps + 15) / abs(math.log10(abs(x)))) + 10
    c = hpl_coeffs(w, N)
    tot = mp.mpf(0)
    xn = x
    for n in range(1, N + 1):
        tot += c[n] * xn
        xn *= x
    return tot

t0 = time.time()
v = H_series((0,1,0,1), mp.mpf('0.77'))
print("series H0101(0.77) =", mp.nstr(v, 30), f"({time.time()-t0:.2f}s)")
t0 = time.time()
v2 = H_series((0,0,1,1), mp.mpf('0.966'))
print("series H0011(0.966) =", mp.nstr(v2, 30), f"({time.time()-t0:.2f}s)")
# sanity vs polylog
print("check H0001:", mp.nstr(H_series((0,0,0,1), mp.mpf('0.77')) - mp.polylog(4, mp.mpf('0.77')), 5))
print("check H00m1=-Li4(-x):", mp.nstr(H_series((0,0,0,-1), mp.mpf('0.77')) + mp.polylog(4, -mp.mpf('0.77')), 5))
print("check H01=Li2:", mp.nstr(H_series((0,1), mp.mpf('0.77')) - mp.polylog(2, mp.mpf('0.77')), 5))
