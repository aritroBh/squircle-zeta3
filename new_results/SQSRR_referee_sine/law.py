# My own transcription of the claimed law C1 (task statement), evaluated exactly with mpmath.
from mpmath import mp, mpf, log, pi
def law(n):
    eps = 1 if n % 4 == 1 else -1
    a, c = n-eps, n+eps
    H = lambda m: sum(mpf(1)/j for j in range(1, m+1))
    tau = pi**2/8 - sum(mpf(1)/(2*i-1)**2 for i in range(1, a//4+1))
    beta = mpf(eps)/((n*n-1)*(2*((n-1)//4)+1))
    rho = -mpf(2*n*(3*n-1))/(a**2*c**3) if eps == 1 else mpf(2*(3*n*n+3*n-2))/(a**2*c**3)
    b2 = (beta*log(2) - mpf(2*eps*n)/(a**2*c)*tau - mpf(2*eps)/a**3*H(a//4)
          + mpf(2*eps*(3*n*n+2*eps*n+3))/(a**3*c**2)*H(a//2) + rho)
    return beta, b2
