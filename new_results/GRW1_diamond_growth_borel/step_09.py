
import mpmath as mp
mp.mp.dps = 60
from fractions import Fraction
from math import factorial
import numpy as np

N = 240
def Hseries_rhombus(q, N):
    a = [(Fraction(1,factorial(n)) if n%2==0 else q/Fraction(factorial(n))) for n in range(N+1)]
    g = [Fraction(0)]*(N+1); g[0]=Fraction(1)
    for k in range(1, N+1):
        g[k] = (a[k] - sum(g[i]*g[k-i] for i in range(1,k)))/2
    H = [Fraction(0)]*(N+1)
    for n in range(N+1):
        H[n] = sum((g[k] + ((-1)**k)*g[k]) * Fraction((-1)**(n-k), 2**(n-k)*factorial(n-k)) for k in range(0, n+1))
    return H

q = Fraction(-3,5)
H = Hseries_rhombus(q, N)
A, B = 1, 2
r = mp.sqrt(5); L = mp.log(2); rho = mp.sqrt(L**2+mp.pi**2/4)
th0 = mp.atan2(mp.pi/2, L)
def bracket(n):
    return B*((-1)**n)*mp.cos((n-mp.mpf('0.5'))*th0) - A*mp.sin((n-mp.mpf('0.5'))*th0)
def lawH(n):
    return -(mp.sqrt(2)/(r*mp.sqrt(mp.pi)))*rho**(mp.mpf('0.5')-n)*n**(-mp.mpf('1.5'))*bracket(n)
Hn = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in H]

safe = [n for n in range(100, 241) if abs(bracket(n)) > 1.0]
print(f"safe n in [100,240] with |bracket|>1: {len(safe)} of 141;  min|bracket| overall: {mp.nstr(min(abs(bracket(n)) for n in range(100,241)),4)}")
ratios = [Hn[n]/lawH(n) for n in safe]
X = np.array([[1.0, 1.0/n, 1.0/n**2, 1.0/n**3] for n in safe])
Y = np.array([float(v) for v in ratios])
coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
resid = np.max(np.abs(X@coef - Y))
print("fit const=%.12f a=%.6f b=%.4f c=%.3f  maxres=%.2e"%(coef[0],coef[1],coef[2],coef[3],resid))
print("worst |ratio-1| over safe n:", mp.nstr(max(abs(v-1) for v in ratios), 5))
print("sample ratios:", [mp.nstr(ratios[i], 12) for i in range(0, len(ratios), 12)])
