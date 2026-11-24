import os
print("CPU count:", os.cpu_count())
# check load
try:
    print("loadavg:", os.getloadavg())
except Exception as e:
    print(e)
# Check sympy exact S(1) for both models at dps 90 (quick verification)
import mpmath as mp
mp.mp.dps = 90
import sympy as sp
S1A = sp.N(2*sp.sqrt(5)/3, 95)
print("\nModel A exact S(1) = 2sqrt5/3 =", S1A)
S1B_expr = -36 + 21*sp.sqrt(3) + sp.Rational(9,13)*sp.sqrt(35515-20504*sp.sqrt(3))
S1B = sp.N(S1B_expr, 95)
print("Model B exact S(1) (sympy expr) =", S1B)
# Compare with transcript claims
h96_B_d0 = mp.mpf('1.075764930490848519873765587084270092790655697826818297768')
trans_B = mp.mpf('1.0757649304908485365928072496899403631687164306641')
print("\nh96 driver d^0 - exact:", mp.nstr(h96_B_d0 - mp.mpf(str(S1B)), 4))
print("transcript 60d - exact:", mp.nstr(trans_B - mp.mpf(str(S1B)), 4))
