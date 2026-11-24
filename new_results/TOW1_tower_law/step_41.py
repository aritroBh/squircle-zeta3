from towlib import *
from genintegrand import S_poly, p_comp_1mt, binom_half
ok = True
for d in range(1, 9):
    for m in range(0, d+1):
        if p_comp_1mt(S_poly(d,m)) != S_poly(d, d-m):
            ok = False; print('S_m(1-t)!=S_{d-m}(t) at', d, m)
print('S_m(1-t) == S_{d-m}(t) for d<=8:', ok)

# Now rho at depths 0,1,2 from the analytic factor A_i(v)
from math import comb
def C(n,k): return comb(n,k) if (0 <= k <= n) else 0
tstar = QSI(Fr(1,2),0,Fr(1,2),0)
def S_tstar(d,m):
    return QSI(Fr(binom_half(d))) * (IQ**d) * Fr((-1)**m * C(d,m))
def Sp_tstar(d,m):
    Bd_prime = QSI(Fr(C(d,m))) * (IQ**d) * Fr((-1)**m) * QSI(2*(d-2*m),0,-2*d,0)
    c = Fr(C(d-1,m)*(d-m-1) - C(d-1,m-1)*(m-1))
    Bdm1 = QSI(c) * (IQ**(d-1)) * Fr((-1)**m)
    return QSI(Fr(binom_half(d))) * Bd_prime + QSI(Fr(binom_half(d-1))) * (2*IQ) * Bdm1

def degS(d, m):
    P = S_poly(d,m)
    return len(P)-1 if P else 0

def A_i_value(d, i):
    """A_i(i) = -(sqrt2/2) 8^d 4^{-D} i^{2d-1-D} Stilde_i(i) (2i)^{2-2d}/((i^2-a^2)(i^2-b^2))."""
    D = degS(d, i)
    St_i = (QSI(4)*IQ)**D * S_tstar(d, i)   # Stilde(i) = (4i)^D S_i(t*)
    denom = (IQ**2 - AQ**2)*(IQ**2 - BQ**2)  # (-1-a^2)(-1-b^2)
    C0 = QSI(0, Fr(-1,2)) * QSI(8**d) / QSI(4**D)   # -(sqrt2/2) 8^d 4^{-D}
    return C0 * (IQ**(2*d-1-D)) * St_i * ((2*IQ)**(2-2*d)) / denom

def rho_closed(key, delta):
    i, j = key; d = i+j
    mult = 2 if i > j else 1
    A = A_i_value(d, i)
    if delta == 0:
        return A.scale(mult)
    if delta == 1:
        return (A * QSI(0,0,-(d-2),0)).scale(mult)   # -i(d-2) A
    if delta == 2:
        r = Sp_tstar(d,i)/S_tstar(d,i)
        # A''/A = (3d-3)/2 - (i/2)(S'/S) - (d-2)^2 ; rho = A''/2
        val = QSI(Fr(3*d-3,2)) - r*(IQ.scale(Fr(1,2))) - QSI((d-2)**2)
        return (A*val).scale(Fr(mult,2))
bad = 0; tested = 0
for d in range(2, 11):
    for key in keys_of_degree(d):
        for delta in (0,1,2):
            tested += 1
            if rho_closed(key, delta) != rho(key, delta):
                bad += 1; print('MISMATCH', key, delta, rho_closed(key,delta), rho(key,delta))
print(f'rho closed forms depths 0,1,2 vs binteg: tested {tested}, bad {bad}')
