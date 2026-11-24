from towlib import *
from genintegrand import S_poly, binom_half
from math import comb
import sympy as sp

# t* = (1+i)/2 in QSI
tstar = QSI(Fr(1,2),0,Fr(1,2),0)
one = QSI(1)
def qpow(q, n): return q**n

def S_at_tstar(d, m):
    """S_m^{(d)}(t*) closed form: binom(1/2,d) i^d (-1)^m binom(d,m)."""
    bh = binom_half(d)
    return QSI(Fr(bh)) * (IQ**d) * Fr((-1)**m * comb(d,m))

def Sp_at_tstar(d, m):
    """S_m'(t*) closed form (derivative w.r.t. t)."""
    i2d = IQ**d
    Bd = Fr(comb(d,m)) * (2**d)
    Bdp = Bd * 2 * Fr(1)  # placeholder
    # B_d'(t*) = binom(d,m) 2^d (i/2)^d (-1)^m * 2[(d-2m) - d i]
    Bd_prime = QSI(Fr(comb(d,m) * 2**d)) * (IQ**d) * Fr((-1)**m) * QSI(2*(d-2*m), 0, -2*d, 0) / QSI(2**d)
    # careful: (i/2)^d = i^d / 2^d; binom*2^d * (i/2)^d = binom * i^d
    Bd_prime = QSI(Fr(comb(d,m))) * (IQ**d) * Fr((-1)**m) * QSI(2*(d-2*m),0,-2*d,0)
    # B_{d-1}(t*) = (i/2)^{d-1} (-1)^m 2^{d-1} [2 binom(d-1,m) - (m-1) binom(d-1,m-1)]
    c = Fr(2*comb(d-1,m) - (m-1)*comb(d-1,m-1))
    Bdm1 = QSI(c) * (IQ**(d-1)) * Fr((-1)**m)   # times (1/2)^{d-1}: fold below
    Bdm1 = Bdm1 / QSI(2**(d-1))
    # S' = binom(1/2,d) B_d' + binom(1/2,d-1) * 2i * B_{d-1}
    return QSI(Fr(binom_half(d))) * Bd_prime + QSI(Fr(binom_half(d-1))) * (2*IQ) * Bdm1

# verify against S_poly evaluated at t* (exact: S_poly has Fr coefficients; t* components in QSI)
def S_poly_eval(d, m, tq):
    P = S_poly(d, m)
    out = ZERO
    for s, c in enumerate(P):
        out = out + QSI(c) * (tq ** s)
    return out
print('S_m(t*) closed form checks:')
ok = True
for d in range(1, 9):
    for m in range(0, d+1):
        a = S_at_tstar(d,m); b = S_poly_eval(d,m,tstar)
        c1 = Sp_at_tstar(d,m)
        # numerical derivative check: (S(t*+h)-S(t*))/h for small exact h? use symmetric exact: deriv via polynomial
        P = S_poly(d,m)
        Pd = [Fr(s)*c for s,c in enumerate(P)][1:]
        c2 = ZERO
        for s, cf in enumerate(Pd):
            c2 = c2 + QSI(cf) * (tstar ** s)
        ok &= (a==b)
        ok &= (c1==c2)
        if a!=b or c1!=c2:
            print('MISMATCH', d, m, a, b, c1, c2)
print('all d<=8:', ok)
