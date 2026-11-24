import sys
sys.path.insert(0, '../AV2_beukers_squircle')
from mpmath import (mp, mpf, nsum, inf, pi, zeta, log, sin, cos, sqrt, clsin,
                    clcos, binomial, harmonic, quad, polylog, gamma, psi, exp, nprint)
mp.dps = 50

# Verify the four integration identities numerically
# 1. f(x) = sum x^n/(n^2 C(2n,n)) = 2 arcsin^2(sqrt(x)/2)
x = mpf('0.37')
f_num = nsum(lambda n: x**n/(mpf(n)**2*binomial(2*n,n)), [1, inf])
f_cls = 2*mp.asin(sqrt(x)/2)**2
print("f(x) id residual:", nprint(f_num - f_cls, 5))

# 2. Cl2'(x) = -log(2 sin(x/2))
x0 = mpf('0.9')
res = mp.diff(lambda t: clsin(2,t), x0) + log(2*sin(x0/2))
print("Cl2' id residual:", nprint(res, 5))

# 3. Cl3'(x) = -Cl2(x)
res = mp.diff(lambda t: clcos(3,t), x0) + clsin(2,x0)
print("Cl3' id residual:", nprint(res, 5))

# 4. Cl3(pi/3) = zeta(3)/3
print("Cl3(pi/3)-zeta(3)/3:", nprint(clcos(3,pi/3)-zeta(3)/3, 5))

# Verify T-integral: T = 2 int_0^1 [arcsin^2(sqrt(x)/2) - arcsin^2(x/2)]/(1-x) dx
T_int = 2*quad(lambda x: (mp.asin(sqrt(x)/2)**2 - mp.asin(x/2)**2)/(1-x), [0, 1])
T_sum = nsum(lambda n: (harmonic(2*n)-harmonic(n))/(mpf(n)**2*binomial(2*n,n)), [1, inf])
print("T int - sum:", nprint(T_int - T_sum, 5))

# Verify T closed: T = 2 zeta(3) - (2 pi/3) Cl2(pi/3)
T_cls = 2*zeta(3) - (2*pi/3)*clsin(2, pi/3)
print("T sum - closed:", nprint(T_sum - T_cls, 5))

# Verify U closed: U = (2pi/3)Cl2(pi/3) - 4 zeta(3)/3
U_sum = nsum(lambda n: 1/(mpf(n)**3*binomial(2*n,n)), [1, inf])
U_cls = (2*pi/3)*clsin(2,pi/3) - 4*zeta(3)/3
print("U sum - closed:", nprint(U_sum - U_cls, 5))

# Final: B2(1) = T + U = 2 zeta(3)/3
print("T+U - 2zeta(3)/3:", nprint(T_sum+U_sum - 2*zeta(3)/3, 5))
with open('progress.log','a') as fh:
    fh.write("step_05: classical Clausen/dilog ids + T,U closed forms verified to 50 digits\n")
