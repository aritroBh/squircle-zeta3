# Diagnose quadrature accuracy for complex p: compare methods/precisions at p = 1+0.09 and 1+0.09i
from mpmath import mp, mpf, mpc, quad, sqrt
mp.dps = 100
p = mpc(1) + mpc('0.09')*mpc(0, 1)
e = 1/p
f = lambda t: sqrt((e*t**(e-1))**2 + (e*(1-t)**(e-1))**2)
r1 = 8*quad(f, [0, mpf(1)/16, mpf(1)/4, mpf(1)/2])
r2 = 8*quad(f, [0, mpf(1)/16, mpf(1)/4, mpf(1)/2], method='tanh-sinh', degree=8)
r3 = 8*quad(f, [0, mpf(1)/2], method='gauss-legendre', degree=6)
print("default :", mp.nstr(r1, 60))
print("ts deg8 :", mp.nstr(r2, 60))
print("gl deg6 :", mp.nstr(r3, 60))
print("diff12:", mp.nstr(abs(r1-r2), 3), " diff13:", mp.nstr(abs(r1-r3), 3), " diff23:", mp.nstr(abs(r2-r3), 3))
