from mpmath import mp, mpf, log, exp, factorial, comb
mp.dps = 40
# C(V) = 2 e^V |log(1-e^{-V})|, claimed decreasing from ~4.2 (V=0.2) to 2
prev = None; mono = True
for i in range(20, 20001, 37):
    V = mpf(i)/1000
    c = 2*exp(V)*abs(log(1-exp(-V)))
    if prev is not None and c > prev: mono = False
    prev = c
print("C(V) monotone decreasing on [0.02,20]:", mono, "; C(0.2) =", mp.nstr(2*exp(mpf('0.2'))*abs(log(1-exp(-mpf('0.2')))),8), "; C(20) ~", mp.nstr(prev,8))
# exact k_j: (1-e^{-V})^{±2 tau} = sum L^j/j! tau^j  => |k_j| = (C(V)e^{-V})^j/j! exactly
print("identity |k_j| = (C(V)e^{-V})^j/j!: exact since k_j = L^j/j!, L = ±2 log(1-e^{-V})")
# combinatorial inequalities used in S3.3
ok1 = all(factorial(n-j)/factorial(n) <= 1/factorial(j) for n in range(1,60) for j in range(1,n+1))
ok2 = all((j+1)**(n-j+1) >= 2**(n-j+1) for n in range(1,60) for j in range(1,n+1))
print("(n-j)!/n! <= 1/j! for all 1<=j<=n<60:", ok1, "; (j+1)^{n-j+1} >= 2^{n-j+1}:", ok2)
# is the ratio (pi i/8) exact? strength_edge / strength_H
from mpmath import pi as mpi, sqrt, j as jj
