import sympy as sp
from sympy import Rational, factorial

# Q(w) = (1-exp(-2w))/(2w) = sum_{m>=0} a_m w^m, a_m = (-2)^m/(m+1)!
K = 44
a = [Rational((-2)**m, factorial(m+1)) for m in range(K+3)]
# b = sqrt(a): b_0 = 1, recursion
b = [Rational(0)]*(K+3)
b[0] = Rational(1)
for m in range(1, K+3):
    b[m] = (a[m] - sum(b[i]*b[m-i] for i in range(1,m)))/2
# E(w) = sum_{k>=0} e_k w^k, e_k = b_{k+2}
e = [b[k+2] for k in range(K+1)]
print("e_0 =", e[0], "= 5/24?", e[0]==Rational(5,24))
print("e_1 =", e[1], "  e_2 =", e[2], "  e_3 =", e[3], "  e_4 =", e[4])
print("float e_1..e_6:", [float(e[k]) for k in range(1,7)])
# psi(s) = |E(is)|^2 = sum c_m s^{2m}, c_m = (-1)^m sum_{i+j=2m} (-1)^j e_i e_j
c = []
for m in range(0, K//2):
    s = sum(((-1)**j)*e[i]*e[2*m-i] for i in range(0, 2*m+1) for j in [2*m-i])
    # careful: sum over i+j=2m of (-1)^j e_i e_j
    s = sum(((-1)**(2*m-i))*e[i]*e[2*m-i] for i in range(0, 2*m+1))
    c.append(((-1)**m)*s)
print("c_0 =", c[0], "= (5/24)^2?", c[0]==Rational(25,576))
print("c_1 =", c[1], "=", float(c[1]))
print("c_2 =", c[2], "=", float(c[2]))
print("c_3 =", c[3], "=", float(c[3]))
print("c_4 =", c[4], "=", float(c[4]))
print("c_5..c_10:", [float(c[m]) for m in range(5,11)])
# sanity: Psi(u) on u=s^2 vs numeric
import mpmath as mp
mp.mp.dps = 40
def Enum(s):
    w = 1j*mp.mpf(s)
    q = (1-mp.e**(-2*w))/(2*w)
    return -(mp.sqrt(q) - 1 + w/2)/mp.mpf(s)**2
for s in [0.1, 0.5, 1.0, 2.0]:
    poly = sum(mp.mpf(float(c[m]))*(mp.mpf(s)**2)**m for m in range(0,20))
    print(s, mp.nstr(poly,20), mp.nstr(abs(Enum(s))**2,20))
with open('progress.log','a') as f:
    f.write(f"step_03: exact E-series e0..e4 = {e[0]},{e[1]},{e[2]},{e[3]},{e[4]}; c1={c[1]}, c2={c[2]}\n")
