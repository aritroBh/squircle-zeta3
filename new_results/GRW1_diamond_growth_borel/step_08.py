
import mpmath as mp
mp.mp.dps = 80
from fractions import Fraction
from math import factorial
import numpy as np

N = 240
def Hseries_rhombus(q, N):
    """H(z;q) = e^{-z/2}[g(z)+g(-z)], g=sqrt(cosh z + q sinh z); exact rational in q."""
    a = [(Fraction(1,factorial(n)) if n%2==0 else q/Fraction(factorial(n))) for n in range(N+1)]
    g = [Fraction(0)]*(N+1); g[0]=Fraction(1)
    for k in range(1, N+1):
        g[k] = (a[k] - sum(g[i]*g[k-i] for i in range(1,k)))/2
    # H = e^{-z/2}(g(z)+g(-z)): even part of g times 2, then times e^{-z/2}
    H = [Fraction(0)]*(N+1)
    for n in range(N+1):
        acc = Fraction(0)
        for k in range(0, n+1):
            gk = g[k] if k%2==0 else -g[k]  # g_k + (-1)^k g_k... careful: g(-z) = sum (-1)^k g_k z^k
        # do directly
        acc = sum((g[k] + ((-1)**k)*g[k]) * Fraction((-1)**(n-k), 2**(n-k)*factorial(n-k)) for k in range(0, n+1))
        H[n] = acc
    return H

q = Fraction(-3,5)   # A=1,B=2
H = Hseries_rhombus(q, N)
A, B = 1, 2
r = mp.sqrt(5); L = mp.log(2); rho = mp.sqrt(L**2+mp.pi**2/4)
th0 = mp.atan2(mp.pi/2, L)
def lawH(n):
    bracket = B*((-1)**n)*mp.cos((n-mp.mpf('0.5'))*th0) - A*mp.sin((n-mp.mpf('0.5'))*th0)
    return -(mp.sqrt(2)/(r*mp.sqrt(mp.pi)))*rho**(mp.mpf('0.5')-n)*n**(-mp.mpf('1.5'))*bracket
Hn = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in H]
print("q=-3/5 (A=1,B=2): rho =", mp.nstr(rho, 15), " theta0 =", mp.nstr(th0, 15))
print("n, H_n/lawH:")
for n in [40, 80, 120, 160, 200, 240]:
    print(f"  n={n}: {mp.nstr(Hn[n]/lawH(n), 20)}")
# fit ratio = 1 + a/n + b/n^2 + c/n^3 + d/n^4 on [120,240]
ns = list(range(120, 241))
X = np.array([[1.0, 1.0/n, 1.0/n**2, 1.0/n**3, 1.0/n**4] for n in ns])
Y = np.array([float(Hn[n]/lawH(n)) for n in ns])
coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
resid = np.max(np.abs(X@coef - Y))
print("fit const=%.12f a=%.6f b=%.4f c=%.3f d=%.2f maxres=%.2e"%(coef[0],coef[1],coef[2],coef[3],coef[4],resid))

# also q=0 square exact 1/n correction: fit per class a_even, a_odd
q0 = Fraction(0)
H0 = Hseries_rhombus(q0, N)
H0n = [mp.mpf(v.numerator)/mp.mpf(v.denominator) for v in H0]
twopi = 2/mp.pi
def lawH0(n):
    return mp.sqrt(2)*mp.cos(mp.pi*(3-2*n)/4)*twopi**n*n**(-mp.mpf('1.5'))
for cls in [0,1]:
    ns = [n for n in range(120, 241) if n%2==cls]
    X = np.array([[1.0, 1.0/n, 1.0/n**2, 1.0/n**3] for n in ns])
    Y = np.array([float(H0n[n]/lawH0(n)) for n in ns])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = np.max(np.abs(X@coef - Y))
    print(f"square class{cls}: const={coef[0]:.12f} a={coef[1]:.8f} b={coef[2]:.5f} c={coef[3]:.4f} res={resid:.1e}")
print("candidates: -3/16 =", -3/16, " 3pi/8 =", 3*float(mp.pi)/8, " pi/8 =", float(mp.pi)/8)
