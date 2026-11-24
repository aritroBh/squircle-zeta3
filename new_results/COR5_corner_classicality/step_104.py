import sys
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from hpl_series import Li31, H_series
mp.mp.dps = 100
x, y = mp.mpf(-1)/2, mp.mpf(-1)
v1 = Li31(x, y)
# independent double series: Li_{3,1}(x,y) = sum_{n>m>=1} x^n y^m/(n^3 m)
def ds(M):
    tot = mp.mpf(0)
    for n in range(2, M):
        for m in range(1, n):
            tot += x**n * y**m / (mp.mpf(n)**3 * m)
    return tot
import math
for M in (200, 400, 800):
    v2 = ds(M)
    # tail estimate: sum_{n>M} |x|^n/(n^3) * (n-1) ~ |x|^M/(M^2(1-|x|))
    print(f"M={M}: ds={mp.nstr(v2, 30)}  |ds-Li31|={mp.nstr(abs(v2-v1),3)}")
print("Li31:", mp.nstr(v1, 40))
# also H_series
print("H_series:", mp.nstr(H_series((0,0,-1,1), mp.mpf(1)/2), 40), " (should equal -Li31)")
