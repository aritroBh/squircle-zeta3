
import sys
sys.path.insert(0, '<project>/research/08_theory/COR8_depth_two_any_angle')
import mpmath as mp
mp.mp.dps = 50
import sympy as sp
from cor8_series import *
c = sp.Symbol('c')
def sLi_n(g, n):
    g0 = g.subs(c, 0)
    dlogg = sdiff(slog_of_rat(g, c))
    if g0 == 0:
        cur = sscale(slog_of_rat(1 - g, c), -1)
        for k in range(2, n+1):
            cur = sintegrate(smul(cur, dlogg))
        return trim(cur)
    if g0 == 1:
        cur = sscale(slog_of_rat(1 - g, c), -1)
        for k in range(2, n+1):
            cur = sintegrate(smul(cur, dlogg))
            cur[(0, 0)] = mp.zeta(k)
        return trim(cur)
    cur = {(0, 0): mp.re(mp.polylog(1, mp.mpf(g0)))}
    for k in range(2, n+1):
        cur = sintegrate(smul(cur, dlogg))
        cur[(0, 0)] = mp.re(mp.polylog(k, mp.mpf(g0)))
    return trim(cur)
g = 1/(1-c)
S2 = sLi_n(g, 2)
print("my Li2(1/(1-c)) series:", sorted([(k, mp.nstr(v, 10)) for (k,l), v in S2.items() if k <= 5]))
# direct: Li2(1/(1-c)) = zeta(2) + int_0^c (-log t + log(1-t))/(1-t) dt
# expand integrand: (-log t)(1+t+t^2+...) + (log(1-t))(1+t+t^2+...)
# -log t * sum t^k: integral -c^{k+1}(log c/(k+1) - 1/(k+1)^2)
# log(1-t)/(1-t) = -(sum t^m/m)(sum t^j) = -sum_n H_n t^n... wait log(1-t)/(1-t) = -sum_{n>=1} H_n t^n? 
# d/dc [1/2 log^2(1-c)] = -log(1-c)/(1-c) = sum H_n c^n
import mpmath as mp2
ref = {}
z2 = mp.zeta(2)
ref[(0,0)] = z2
H = 0
for n in range(1, 12):
    H += mp.mpf(1)/n   # H_n
# integrand series: A(t) = -log t * sum_{k>=0} t^k + (-sum_{n>=1} H_n t^n)?? check: -log(1-t)/(1-t) = sum_{n>=0} H_{n+1}... 
# 1/2 log^2(1-t) = sum_{n>=2} H_{n-1}/n t^n; derivative: -log(1-t)/(1-t) = sum_{n>=2} H_{n-1} t^{n-1} = sum_{m>=1} H_m t^m
A = {}
for k in range(0, 10):
    A[(k,1)] = A.get((k,1),0) - 1      # -log t * t^k
    if k >= 1:
        Hk = sum(mp.mpf(1)/m for m in range(1, k+1))
        A[(k,0)] = A.get((k,0),0) - Hk  # -log(1-t)/(1-t) = sum H_k t^k
I = sintegrate(A)
I[(0,0)] = z2
print("reference:", sorted([(k, mp.nstr(v, 10)) for (k,l), v in I.items() if k <= 5]))
