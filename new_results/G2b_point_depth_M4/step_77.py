
import mpmath as mp
mp.mp.dps = 60
from fractions import Fraction as Fr

def prepend_a(s, a, N):
    t = [mp.mpf(0)] * (N + 1)
    aP = mp.mpf(1); acc = mp.mpf(0)
    for n in range(1, N):
        aP *= a; acc += s[n] * aP
        t[n + 1] = -acc / (aP * a * (n + 1))
    return t

def J_series(word, N=2000):
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

print("identity J-word == +Li32(x,y) (paper), six convergent atoms:")
for (x, y) in [(Fr(1,3), Fr(1)), (Fr(2,3), Fr(1)), (Fr(1,2), Fr(2,3)), (Fr(1,3), Fr(3,2)), (Fr(8,9), Fr(3,4)), (Fr(-1,2), Fr(1))]:
    jv = J_series((0, 0, mp.mpf(1)/x, 0, mp.mpf(1)/(x*y)))
    lv = li32(x.numerator/x.denominator, y.numerator/y.denominator)
    print(f"  x={str(x):>4}, y={str(y):>3}: |J - Li32| = {mp.nstr(abs(jv-lv),3)}")
