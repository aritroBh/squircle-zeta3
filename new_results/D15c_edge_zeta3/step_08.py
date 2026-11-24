
import mpmath as mp
print("dps before:", mp.mp.dps)
mp.mp.dps = 60
import sympy as sp
print("dps after sympy import:", mp.mp.dps)

# Structural verification: D3 is a polynomial in the logs of total degree <=3,
# because S0 = T0 = sum(s_i) is log-free and all denominators are powers of S0 and B=1+Q0^2 (log-free).
# => G21 extracted by step_09's degree-3 fit equals the TRUE coefficient c21 of l0^2*l1. 
# Let me confirm the max log-degree of each piece symbolically at a numeric level first (m=4).
from poly_coeffs import make_coeffs, regular_zeros
mp.mp.dps = 60
print("dps after poly_coeffs import:", mp.mp.dps)
m=4
zeros = regular_zeros(m)
coeffs = make_coeffs(zeros,[1]*m,[0]*m)
th = mp.mpf('0.3')
print("D3 sample (m=4, th=0.3):", mp.nstr(coeffs(th)[3],20))
