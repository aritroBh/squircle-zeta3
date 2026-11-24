import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 80
from oct5d3_atoms import li32_xy
r2=mp.sqrt(2); b=r2-1
def ds(x,y,N):
    tot = mp.mpf(0); Hy = mp.mpf(0); yp=mp.mpf(1)
    for n in range(1,N+1):
        yp *= y; Hy += yp/(n*n)
        if n>=2: tot += x**n/mp.mpf(n)**3*Hy
    return tot
v1 = li32_xy(b,-1)
v2 = ds(b,-1,400)
v3 = ds(b,-1,800)
print("li32_xy:", mp.nstr(v1, 40))
print("ds N=400:", mp.nstr(v2, 40))
print("ds N=800:", mp.nstr(v3, 40))
print("v1-v2:", mp.nstr(v1-v2, 5), " v2-v3:", mp.nstr(v2-v3, 5))
