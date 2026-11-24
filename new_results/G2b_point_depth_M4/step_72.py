
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 60

def J_li32_word(x, y, N=2500):
    """J((0,0,1/x,0,1/(xy)); 0->1) by power series in t (letters 1/x, 1/(xy) with abs >= 4/3)."""
    a, b = mp.mpf(1)/x, mp.mpf(1)/(x*y)
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = -(mp.mpf(1)/b)**m / m
    for letter in (0, a, 0, 0):
        t = [mp.mpf(0)] * (N + 1)
        if letter == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N):
                acc = acc / letter + s[m]
                t[m + 1] = acc / (m + 1)
        s = t
    return mp.fsum(s[1:])

def li32_paper(x, y, N=2500):
    x = mp.mpf(x); y = mp.mpf(y)
    S = mp.mpf(0); ym = y; xn = x; tot = mp.mpf(0)
    for n in range(2, N + 1):
        S += ym / ((n-1)*(n-1)); ym *= y; xn *= x
        tot += xn * S / (n*n*n)
    return tot

print("off-by-x check: J((0,0,1/x,0,1/(xy));0->1) =?= -Li_{3,2}(x,y) (paper convention)")
for (x, y) in [(Fr(1,3), Fr(1)), (Fr(2,3), Fr(1)), (Fr(1,2), Fr(2,3)), (Fr(1,3), Fr(3,2)), (Fr(8,9), Fr(3,4)), (Fr(-1,2), Fr(1))]:
    jv = J_li32_word(x, y)
    lv = li32_paper(x.numerator/x.denominator, y.numerator/y.denominator)
    print(f"  x={str(x):>4}, y={str(y):>3}: diff |J + Li32| = {mp.nstr(abs(jv+lv),3)}")
