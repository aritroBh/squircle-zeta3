import os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
from mpmath import mp, mpf, cos, pi, log, zeta, harmonic, quad
mp.dps = 50
from math import comb

def catalan(m): return comb(2*m, m)//(m+1)
def Kjv(j, x):
    return sum((-2*catalan(m-1))*(-1)**(j-m)*comb(m+j-2, j-m)*x**m for m in range(1, j+1))

def Hgen(r, m):  # H^{(r)}_m
    return sum(mpf(1)/i**r for i in range(1, m+1))

def Lam(m, j):
    S = j**m + (-1)**m*(j-1)**m - 1
    return (-1)**m*zeta(m)*S/m + Hgen(m, j-1)*(j**m-(j-1)**m)/m

def bellY(r, j):  # [u^r] exp(sum_{m>=1} Lam_m u^m)
    # complete Bell polynomial Y_r(Lam_1..Lam_r)/r!
    Ls = [Lam(m, j) for m in range(1, r+1)]
    # recurrence: Y_0=1; Y_r = (1/r) sum_{i=1}^r i*Ls[i-1]*Y_{r-i}  ... for exp sum: y_r = sum_{i=1}^r (i/r) L_i y_{r-i}
    Y = [mpf(1)]
    for rr in range(1, r+1):
        Y.append(sum(Ls[i-1]*Y[rr-i]*i for i in range(1, rr+1))/rr)
    return Y[r]

def Mk_formula(k, alpha, J=600):
    c2 = cos(alpha/2)**2
    return sum(Kjv(j, c2)*(-1)**(j-1)/j*bellY(k+1, j) for j in range(1, J+1))

# test M_0
for alpha in [mpf('2.5'), mpf('2.0'), mpf('1.2')]:
    c = cos(alpha/2)
    M0 = -(1+c)*log(1+c)-(1-c)*log(1-c)
    f0 = Mk_formula(0, alpha)
    print('M0', mp.nstr(alpha,6), mp.nstr(M0, 15), mp.nstr(f0, 15), 'diff', mp.nstr(M0-f0, 5))
