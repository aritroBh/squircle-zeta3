import os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
from mpmath import mp, mpf, cos, sin, pi, gamma, log, sqrt, quad, binomial as bnm
mp.dps = 50

def M_int(u, alpha, J=2000):
    # original corner integral (paper): int_0^1 [(1+g)^(-1-u) sqrt(1+g^(2-2u) - 2 cos a g^(1-u)) - 1] g^(u-1) dg
    ca = cos(alpha)
    f = lambda g: ((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*ca*g**(1-u)) - 1)*g**(u-1)
    # split for the -1 tail: int_0^1 [..] dg with g^{u-1}; singularity at 0 integrable
    return quad(f, [0, 1])

def B_j(u, j):
    v = 1-u
    return gamma(1+(j-1)*v)*gamma(1-j*v)/gamma(1+u)

def catalan(m):
    from math import comb
    return comb(2*m, m)//(m+1)

def M_betasum(u, alpha, J=400):
    c2 = cos(alpha/2)**2
    v = 1-u
    tot = mpf(0)
    for j in range(1, J+1):
        Kj = sum((-2*catalan(n-1))*(-1)**(j-n)*bnm(n+j-2, 2*n-2)*c2**n for n in range(1, j+1))
        tot += Kj*B_j(u, j)
    return tot

for alpha in [mpf('2.5'), mpf('2.9')]:
    for u in [mpf('0.2'), mpf('0.35')]:
        a = M_int(u, alpha)
        b = M_betasum(u, alpha)
        print(f'alpha={alpha}, u={u}: int={mp.nstr(a,12)} betasum={mp.nstr(b,12)} ratio={mp.nstr(a/b,8)}')
