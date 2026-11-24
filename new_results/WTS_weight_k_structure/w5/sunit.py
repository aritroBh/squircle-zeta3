# enumerate z = +-(1+sqrt2)^m sqrt2^n with 1-z also +-(1+sqrt2)^m' sqrt2^n' (S-unit equation in Q(sqrt2), S={sqrt2,infty}); exact Fractions
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as F
def mul(x, y): return (x[0]*y[0] + 2*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def pw(x, n):
    r = (F(1), F(0))
    if n < 0:
        # inverse of (p,q): (p - q sqrt2)/(p^2 - 2q^2)
        N = x[0]**2 - 2*x[1]**2; x = (x[0]/N, -x[1]/N); n = -n
    for _ in range(n): r = mul(r, x)
    return r
a = (F(1), F(1)); s2 = (F(0), F(1))
def is_sunit(x):
    if x == (0, 0): return False
    # N = p^2 - 2 q^2 must be +-2^k
    N = x[0]**2 - 2*x[1]**2
    if N < 0: N = -N
    num, den = N.numerator, N.denominator
    for v in (num, den):
        while v % 2 == 0: v //= 2
        if v != 1: return False
    return True
sols = []
R = 8
for m in range(-R, R+1):
    for n in range(-R, R+1):
        for sg in (1, -1):
            z = pw(a, m); z = mul(z, pw(s2, n)); z = (sg*z[0], sg*z[1])
            w = (1 - z[0], -z[1])
            if is_sunit(w): sols.append((m, n, sg, z, w))
import math
for m, n, sg, z, w in sols:
    val = float(z[0]) + float(z[1])*math.sqrt(2)
    print(f"z={'-' if sg<0 else '+'}a^{m} s2^{n} = {float(z[0]):.4g}+{float(z[1]):.4g}sqrt2 = {val:.6f}   1-z={float(w[0]):.4g}+{float(w[1]):.4g}sqrt2")
