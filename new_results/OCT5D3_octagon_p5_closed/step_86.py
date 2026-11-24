import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 80
from oct5d3_atoms import li32_xy, k32_value
r2=mp.sqrt(2); b=r2-1
def ds_correct(x,y,N):
    tot = mp.mpf(0); Hy = mp.mpf(0); yp=mp.mpf(1); xn = x*x
    for n in range(1,N+1):
        yp *= y; Hy += yp/(n*n)
        if n+1>=2 and n+1<=N: pass
    # rebuild properly
    tot = mp.mpf(0); Hy = mp.mpf(0); yp=mp.mpf(1)
    for m in range(1, N):
        yp *= y; Hy += yp/(m*m)
        n = m+1
        tot += x**n/mp.mpf(n)**3*Hy
    return tot
for (x,y) in [(b,-1),(b,1),(-b,-1),(-b,1),(-b*b,1)]:
    print(x,y, mp.nstr(abs(li32_xy(x,y)-ds_correct(x,y,1200)),3))
print("K32 =", mp.nstr(k32_value(), 30))
