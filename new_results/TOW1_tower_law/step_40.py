from towlib import *
from genintegrand import S_poly, binom_half
from math import comb
def C(n,k): return comb(n,k) if (0 <= k <= n) else 0
tstar = QSI(Fr(1,2),0,Fr(1,2),0)
def S_at_tstar(d, m):
    return QSI(Fr(binom_half(d))) * (IQ**d) * Fr((-1)**m * C(d,m))
def Sp_at_tstar(d, m):
    Bd_prime = QSI(Fr(C(d,m))) * (IQ**d) * Fr((-1)**m) * QSI(2*(d-2*m),0,-2*d,0)
    c = Fr(C(d-1,m)*(d-m-1) - C(d-1,m-1)*(m-1))
    Bdm1 = QSI(c) * (IQ**(d-1)) * Fr((-1)**m)   # (i/2)^{d-1} 2^{d-1} = i^{d-1}
    return QSI(Fr(binom_half(d))) * Bd_prime + QSI(Fr(binom_half(d-1))) * (2*IQ) * Bdm1
ok = True
for d in range(1, 11):
    for m in range(0, d+1):
        P = S_poly(d,m)
        b = ZERO
        for s, cf in enumerate(P): b = b + QSI(cf)*(tstar**s)
        Pd = [Fr(s)*cf for s,cf in enumerate(P)][1:]
        c2 = ZERO
        for s, cf in enumerate(Pd): c2 = c2 + QSI(cf)*(tstar**s)
        if not (S_at_tstar(d,m)==b and Sp_at_tstar(d,m)==c2):
            ok=False; print('MISMATCH',d,m, Sp_at_tstar(d,m), c2)
print('S(t*), Sp(t*) closed forms verified d<=10:', ok)
