
import mpmath as mp
mp.mp.dps = 80

def Jval0(w, x, N=800):
    w = list(w)
    nz = sum(1 for a in w if a == 1)
    tz = 0
    while w and w[-1] == 0:
        tz += 1; w.pop()
    if not w:
        return mp.log(x) ** tz / mp.factorial(tz)   # all-zero word: log^k/k! (sign in nz=0)
    a = w[-1]
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = mp.mpf(1) / m if a == 1 else mp.mpf((-1) ** (m + 1)) / m
    for a in reversed(w[:-1]):
        t = [mp.mpf(0)] * (N + 1)
        if a == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N):
                acc = acc * a + s[m]; t[m + 1] = acc / (m + 1)
        s = t
    for _ in range(tz):
        s = [mp.mpf(0)] + [s[m] / m for m in range(1, N)] + [mp.mpf(0)]
    Hv = mp.fsum(s[m] * x ** m for m in range(1, N + 1))
    return ((-1) ** nz) * Hv

def hpl_series_val(w, x, N):
    a = w[-1]
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = mp.mpf(1) / m if a == 1 else mp.mpf((-1) ** (m + 1)) / m
    for a in reversed(w[:-1]):
        t = [mp.mpf(0)] * (N + 1)
        if a == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N):
                acc = acc * a + s[m]; t[m + 1] = acc / (m + 1)
        s = t
    return mp.fsum(s[m] * x ** m for m in range(1, N + 1))

print("J((0,0,0);1/2) =", mp.nstr(Jval0((0,0,0), mp.mpf(1)/2), 15), "expected", mp.nstr(mp.log(mp.mpf(1)/2)**3/6, 15))
print("J((0,0,1);1/2) =", mp.nstr(Jval0((0,0,1), mp.mpf(1)/2), 15), "vs -hpl:", mp.nstr(-hpl_series_val((0,0,1), mp.mpf(1)/2, 800), 15))
print("J((1,0,0);1/2) =", mp.nstr(Jval0((1,0,0), mp.mpf(1)/2), 15), "vs -hpl:", mp.nstr(-hpl_series_val((1,0,0), mp.mpf(1)/2, 800), 15))
print("J((0,1,0);1/2) =", mp.nstr(Jval0((0,1,0), mp.mpf(1)/2), 15), "vs -hpl:", mp.nstr(-hpl_series_val((0,1,0), mp.mpf(1)/2, 800), 15))
