
import mpmath as mp
from fractions import Fraction as Fr
import pickle, math, time
mp.mp.dps = 80
from sympy.functions.combinatorial.numbers import stirling

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

def I_int(m, j):
    return mp.quad(lambda s: s**m * mp.log(1+mp.e**(-s))**j, [0, mp.inf])

def I_stir(m, j, L=3000):
    tot = mp.mpf(0)
    for l in range(j, L+1):
        s = int(stirling(l, j))  # signed
        tot += mp.mpf(s) / (mp.factorial(l) * mp.mpf(l)**(m+1))
    return mp.factorial(j)*mp.factorial(m)*tot

def I_euler(m, j, N=4000):
    tot = mp.mpf(0)
    e = [mp.mpf(1)] + [mp.mpf(0)]*(j-1)
    for n in range(1, N+1):
        if n >= j:
            tot += ((-1)**n) * e[j-1] / mp.mpf(n)**(m+2)
        for k in range(min(j-1, n), 0, -1):
            e[k] += e[k-1]/n
    return mp.factorial(m) * ((-1)**j) * tot

print("m j  integral / stirling / euler, differences")
for (m,j) in [(0,1),(1,1),(2,1),(0,2),(1,2),(3,1),(2,2),(0,3),(4,2),(0,4),(6,2),(0,6),(7,1),(0,8)]:
    a = I_int(m,j); b = I_stir(m,j); c = I_euler(m,j)
    print(m, j, mp.nstr(a,22))
    print("   |a-b|", mp.nstr(abs(a-b),4), " |a-c|", mp.nstr(abs(a-c),4))
