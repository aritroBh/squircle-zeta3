
import mpmath as mp
mp.mp.dps = 60

def prepend_a(s, a, N):
    t = [mp.mpf(0)] * (N + 1)
    aP = mp.mpf(1); acc = mp.mpf(0)
    for n in range(1, N):
        aP *= a; acc += s[n] * aP
        t[n + 1] = -acc / (aP * a * (n + 1))
    return t

def J_series(word, N=2000):
    """J(word; 0->1), letters nonzero or 0, rightmost nonzero, all |nonzero letters| >= 1."""
    w = list(word)
    a = w[-1]
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = -(mp.mpf(1)/a)**m / m
    for letter in reversed(w[:-1]):
        if letter == 0:
            s = [mp.mpf(0)] + [s[m]/m for m in range(1, N)] + [mp.mpf(0)]
        else:
            s = prepend_a(s, mp.mpf(letter), N)
    return mp.fsum(s[1:])

def li32(x, y, N=2000):
    x = mp.mpf(x); y = mp.mpf(y)
    S = mp.mpf(0); ym = y; xn = x; tot = mp.mpf(0)
    for n in range(2, N + 1):
        S += ym / ((n-1)*(n-1)); ym *= y; xn *= x
        tot += xn * S / (n*n*n)
    return tot

J = J_series((0, 0, 2, 0, 3))
print("J((0,0,2,0,3);1) series =", mp.nstr(J, 22))
v = li32(mp.mpf(1)/2, mp.mpf(2)/3)
print(" Li32(1/2,2/3)      =", mp.nstr(v, 22))
print("diff J+Li32 =", mp.nstr(abs(J + v), 3), "  diff J-Li32 =", mp.nstr(abs(J - v), 3))
