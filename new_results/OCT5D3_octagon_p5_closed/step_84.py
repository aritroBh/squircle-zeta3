import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
src = open("oct5d3_atoms.py").read()
addition = '''

def li32_xy(x, y):
    """Li_{3,2}(x,y) = sum_{n>m>0} x^n y^m/(n^3 m^2), |x|<1, |y|<=1."""
    tot = mp.mpf(0); Hy = mp.mpf(0); yp = mp.mpf(1); xp = x; n = 1
    thresh = mp.mpf(10)**(-mp.mp.dps-15)
    while True:
        yp *= y; Hy += yp/(n*n)
        n += 1
        xp *= x
        tot += xp*Hy/mp.mpf(n)**3
        if n > 10 and abs(xp) < thresh: break
    return tot

def k32_value():
    r2 = mp.sqrt(2); b = r2-1
    return (li32_xy(-b*b, 1) - 8*li32_xy(b, 1) - 32*li32_xy(-b, 1)
            - 16*li32_xy(b, -1) - 16*li32_xy(-b, -1))
'''
open("oct5d3_atoms.py","a").write(addition)
import subprocess, sys
r = subprocess.run([sys.executable, "-c", """
import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 100
from oct5d3_atoms import k32_value, li32_xy
print("K32 value =", mp.nstr(k32_value(), 40))
# sanity: Li32(x,-1) vs direct double series
def ds(x,y,N=2000):
    tot = mp.mpf(0)
    for n in range(2,N+1):
        tot += x**n/n**3*sum(y**m/m**2 for m in range(1,n))
    return tot
r2=mp.sqrt(2); b=r2-1
print("check li32_xy(b,-1):", mp.nstr(abs(li32_xy(b,-1)-ds(b,-1)),3))
"""], capture_output=True, text=True)
print(r.stdout, r.stderr[-500:] if r.returncode else "")
