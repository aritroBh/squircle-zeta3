from towlib import *
from genintegrand import S_poly, binom_half
from math import comb
def C(n,k):
    return comb(n,k) if (0 <= k <= n) else 0
tstar = QSI(Fr(1,2),0,Fr(1,2),0)
def S_at_tstar(d, m):
    return QSI(Fr(binom_half(d))) * (IQ**d) * Fr((-1)**m * C(d,m))
def Sp_at_tstar(d, m):
    Bd_prime = QSI(Fr(C(d,m))) * (IQ**d) * Fr((-1)**m) * QSI(2*(d-2*m),0,-2*d,0)
    c = Fr(2*C(d-1,m) - (m-1)*C(d-1,m-1))
    Bdm1 = QSI(c) * (IQ**(d-1)) * Fr((-1)**m) / QSI(2**(d-1))
    return QSI(Fr(binom_half(d))) * Bd_prime + QSI(Fr(binom_half(d-1))) * (2*IQ) * Bdm1
ok = True
for d in range(1, 10):
    for m in range(0, d+1):
        P = S_poly(d,m)
        a = S_at_tstar(d,m)
        b = ZERO
        for s, cf in enumerate(P): b = b + QSI(cf)*(tstar**s)
        c1 = Sp_at_tstar(d,m)
        Pd = [Fr(s)*cf for s,cf in enumerate(P)][1:]
        c2 = ZERO
        for s, cf in enumerate(Pd): c2 = c2 + QSI(cf)*(tstar**s)
        if not (a==b and c1==c2):
            ok = False; print('MISMATCH', d, m)
print('S(t*), S'(t*') closed forms verified d<=9:', ok)
