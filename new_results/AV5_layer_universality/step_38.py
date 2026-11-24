import mpmath
from mpmath import mp, mpf, log, exp, quad, zeta, nstr
mp.dps = 130
A = quad(lambda x: log(1+x)**2/x, [0, mpf('0.5'), 1])
print("int_0^1 log(1+x)^2/x =", nstr(A, 60))
print("zeta(3)/4            =", nstr(zeta(3)/4, 60))
print("diff:", nstr(A - zeta(3)/4, 6))
B = quad(lambda s: log(1+exp(-2*s))**2, [0, 5, 20, 40])
print("I(0,2) =", nstr(B,60))
print("A/2    =", nstr(A/2,60))
print("diff:", nstr(B-A/2,6))
