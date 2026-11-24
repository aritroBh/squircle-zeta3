
import mpmath as mp
from fractions import Fraction as Fr
import pickle, math, time
mp.mp.dps = 80

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

# I_{m,j} = int_0^inf s^m log^j(1+e^{-s}) ds   (definition)
def I_int(m, j):
    return mp.quad(lambda s: s**m * mp.log(1+mp.e**(-s))**j, [0, mp.inf])

# Stirling series:  I_{m,j} = j! m! sum_{l>=j} s_signed(l,j)/(l! l^{m+1})
from sympy import stirling
def I_stir(m, j, L=3000):
    tot = mp.mpf(0)
    for l in range(j, L+1):
        s = int(stirling(l, j))  # signed
        tot += mp.mpf(s) / (mp.factorial(l) * mp.mpf(l)**(m+1))
    # tail bound: |s(l,j)|/l! <= (log l)^{j-1}/l roughly; terms alternate-decrease after a while; add remainder estimate
    return mp.factorial(j)*mp.factorial(m)*tot

# Euler sum:  I_{m,j} = m! * (-1)^j * zeta(bar(m+2), 1^{j-1}),  zeta(bar s,1..) = sum_{n1>...>nj} (-1)^{n1}/(n1^s n2...nj)
def I_euler(m, j, N=4000):
    # sum_{n} (-1)^n H^{(j-1)}_{n-1} / n^{m+2}  where H^{(j-1)}_{n-1} = |s(n,j)|/(n-1)!
    # compute |s(n,j)|/(n-1)! = e_{j-1}(1,1/2,...,1/(n-1))  recursively
    tot = mp.mpf(0); sgn = 1
    e = [mp.mpf(1)] + [mp.mpf(0)]*(j-1)   # e[k] = elem sym of first n-1 reciprocals
    for n in range(1, N+1):
        if n >= j:
            tot += ((-1)**n) * e[j-1] / mp.mpf(n)**(m+2)
        # update e with 1/n
        for k in range(min(j-1, n), 0, -1):
            e[k] += e[k-1]/n
    return mp.factorial(m) * ((-1)**j) * tot

print("m j  integral                stirling(L=3000)        euler(N=4000)")
for (m,j) in [(0,1),(1,1),(2,1),(0,2),(1,2),(3,1),(2,2),(0,3),(4,2),(0,4),(6,2),(0,6),(7,1),(0,8)]:
    a = I_int(m,j); b = I_stir(m,j); c = I_euler(m,j)
    print(m, j, mp.nstr(a,25), mp.nstr(b,25), mp.nstr(c,25))
    print("   |a-b|", mp.nstr(abs(a-b),5), " |a-c|", mp.nstr(abs(a-c),5))
