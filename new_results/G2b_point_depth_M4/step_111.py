
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 80
from g2b_dbg import s05
C0 = Fr(1,2)

def hpl_series_val(w, x, N):
    """H_w(x) for w ending NONZERO; letters 0,+-1."""
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

def Jval0(w, x, N=800):
    """J(w; 0->x), letters {0,+-1}, regularized; = (-1)^#1 H_w(x)."""
    w = list(w)
    nz = sum(1 for a in w if a == 1)
    # strip trailing zeros (they correspond to outer 1/t integrations)
    tz = 0
    while w and w[-1] == 0:
        tz += 1; w.pop()
    if not w:
        Hv = mp.log(x) ** tz / mp.factorial(tz)
    else:
        s = None
        # build series of H_w(x)
        a = w[-1]
        ss = [mp.mpf(0)] * (N + 1)
        for m in range(1, N + 1):
            ss[m] = mp.mpf(1) / m if a == 1 else mp.mpf((-1) ** (m + 1)) / m
        for a in reversed(w[:-1]):
            t = [mp.mpf(0)] * (N + 1)
            if a == 0:
                for m in range(1, N + 1): t[m] = s[m] / m if s else 0
            else:
                acc = mp.mpf(0)
                for m in range(0, N):
                    acc = acc * a + s[m]; t[m + 1] = acc / (m + 1)
            s = t
        # trailing zeros: divide coefficients by m
        for _ in range(tz):
            s = [mp.mpf(0)] + [s[m] / m for m in range(1, N)] + [mp.mpf(0)]
        Hv = mp.fsum(s[m] * x ** m for m in range(1, N + 1))
    return ((-1) ** nz) * Hv

# sanity: J((0,0,0);0->1/2) = log^3(1/2)/6
print("J((0,0,0);1/2) =", mp.nstr(Jval0((0,0,0), mp.mpf(1)/2), 15), " expected", mp.nstr(mp.log(mp.mpf(1)/2)**3/6, 15))
print("J((0,0,1);1/2) =", mp.nstr(Jval0((0,0,1), mp.mpf(1)/2), 15), "vs hpl", mp.nstr(-hpl_series_val((0,0,1), mp.mpf(1)/2, 800), 15))
print("J((0,1,0);1/2) =", mp.nstr(Jval0((0,1,0), mp.mpf(1)/2), 15))

# delta_num for the Li5 word with fixed Jval0
exec(open("/dev/stdin").read()) if False else None
