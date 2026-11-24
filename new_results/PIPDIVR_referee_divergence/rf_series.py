# Truncated power-series helpers on lists of flint arb/acb (own code for the PIPDIVR referee).
from flint import arb, acb

def mul(a, b, N):
    out = [0 * a[0] for _ in range(N)]
    for i, ai in enumerate(a[:N]):
        if ai == 0:
            continue
        for j in range(min(len(b), N - i)):
            out[i + j] += ai * b[j]
    return out

def inv(a, N):
    out = [1 / a[0]]
    for n in range(1, N):
        s = 0 * a[0]
        for k in range(1, min(n, len(a) - 1) + 1):
            s += a[k] * out[n - k]
        out.append(-s / a[0])
    return out

def deriv(a):
    return [a[i] * i for i in range(1, len(a))]

def integ(a, c0):
    return [c0] + [a[i] / (i + 1) for i in range(len(a))]

def log(a, N, c0):
    # log of series a with a[0] != 0; c0 = chosen value of log(a[0])
    d = mul(deriv(a) + [0 * a[0]], inv(a, N), N - 1)
    return integ(d, c0)[:N]

def exp(a, N):
    # exp of series a, via E' = a' E
    e0 = a[0].exp()
    da = deriv(a) + [0 * a[0]]
    out = [e0]
    for n in range(1, N):
        s = 0 * a[0]
        for k in range(1, n + 1):
            if k - 1 < len(da):
                s += da[k - 1] * out[n - k]
        out.append(s / n)
    return out
