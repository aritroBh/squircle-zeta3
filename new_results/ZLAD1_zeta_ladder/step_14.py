from mpmath import mp, mpf, cos, pi, log, zeta, quad
mp.dps = 50
from math import comb
def catalan(m): return comb(2*m, m)//(m+1)
def Kjv(j, x):
    return sum((-2*catalan(m-1))*(-1)**(j-m)*comb(m+j-2, j-m)*x**m for m in range(1, j+1))
def Hgen(r, m):
    return sum(mpf(1)/i**r for i in range(1, m+1))
def Lam(m, j):
    S = j**m + (-1)**m*(j-1)**m - 1
    zt = (-1)**m*zeta(m)*S/m if S != 0 else mpf(0)
    return zt + Hgen(m, j-1)*(j**m-(j-1)**m)/m
def bellY(r, j):
    Ls = [Lam(m, j) for m in range(1, r+1)]
    Y = [mpf(1)]
    for rr in range(1, r+1):
        Y.append(sum(Ls[i-1]*Y[rr-i]*i for i in range(1, rr+1))/rr)
    return Y[r]
def Mk_formula(k, alpha, J=800):
    c2 = cos(alpha/2)**2
    return sum(Kjv(j, c2)*(-1)**(j-1)/j*bellY(k+1, j) for j in range(1, J+1))

for alpha in [mpf('2.5'), mpf('2.0'), mpf('1.2')]:
    c = cos(alpha/2)
    M0 = -(1+c)*log(1+c)-(1-c)*log(1-c)
    f0 = Mk_formula(0, alpha)
    print('M0', mp.nstr(alpha,6), mp.nstr(M0, 12), mp.nstr(f0, 12), 'diff', mp.nstr(M0-f0, 4))

# and M_k for k>=1 via the corner integral
def M_int_k(k, alpha):
    ca = cos(alpha)
    import mpmath as M
    u = M.mpf('0.0')
    # M(u;alpha) = sum M_k u^k; compute series coefficient numerically via Cauchy? simpler: finite differences unstable.
    # Use the integral representation with small-u series of integrand? Instead: numeric series extraction.
    return None
# quick check of M_k at alpha=pi/2 against known squircle M_k (paper): M_0(pi/2)=?
alpha = pi/2
c = cos(alpha/2)
print('M0(pi/2) =', mp.nstr(-(1+c)*log(1+c)-(1-c)*log(1-c), 15), ' formula:', mp.nstr(Mk_formula(0, alpha), 15))
