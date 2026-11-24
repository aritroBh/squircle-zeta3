"""w6 identification step 26: two small identities used in the write-up (all printed numbers from this script).
 (i) Taylor coefficients of log E(y), E(y)=2^y Gamma(1-y) Gamma(1+y/2)/Gamma(1-y/2), equal s_n/n with s_n = eta(n) (n odd), zeta(n) (n even);
 (ii) series of (1 - sinh z)/cosh(z)^2 equals sum_m eps_m A_{m+1} z^m / m!  (A_n Euler zigzag numbers, eps_m = +,-,-,+ periodic)."""
import sys
sys.dont_write_bytecode = True
import mpmath as mp
import sympy as sp
mp.mp.dps = 60
print("mp.mp.dps =", mp.mp.dps)
logE = lambda y: y * mp.log(2) + mp.loggamma(1 - y) + mp.loggamma(1 + y / 2) - mp.loggamma(1 - y / 2)
tc = mp.taylor(logE, 0, 10)
for n in range(1, 11):
    s = (1 - mp.power(2, 1 - n)) * mp.zeta(n) if n % 2 == 1 and n > 1 else (mp.log(2) if n == 1 else mp.zeta(n))
    print(f"n={n}: [y^n] log E - s_n/n = {mp.nstr(tc[n] - s / n, 3)}")
z = sp.symbols("z")
ser = sp.series((1 - sp.sinh(z)) / sp.cosh(z) ** 2, z, 0, 10).removeO()
x = sp.symbols("x")
zig = sp.series(sp.sec(x) + sp.tan(x), x, 0, 12).removeO()
A = [int(zig.coeff(x, n) * sp.factorial(n)) for n in range(11)]
eps = [1, -1, -1, 1]
ok = all(sp.simplify(ser.coeff(z, m) - sp.Rational(eps[m % 4] * A[m + 1], sp.factorial(m))) == 0 for m in range(10))
print("zigzag numbers A_0..A_10 =", A)
print("series of (1-sinh z)/cosh^2 z equals sum eps_m A_{m+1} z^m/m! for m=0..9:", ok)
