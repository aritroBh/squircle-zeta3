from towlib import *
from math import comb
def Gm_ipos(i,j,k):
    """+i principal part coefficients e_m (m>=1) of Gm_ij^(k)."""
    H = tower(k)[(i,j)]
    maxr = max(H) if H else 1
    out = {}
    for r in range(2, maxr+1):
        e = unq4(H.get(r,(0,0,0,0)))
        out[r-1] = e.scale(Fr(-1, r-1))
    return out

def Psi_orders(i,j,k):
    """orders of +i principal part of Gm(v) - Gm(1/v): psi_r = e_r - sum_{m>=r} conj(e_m)(-i)^m binom(m,r) i^{m-r}."""
    e = Gm_ipos(i,j,k)
    out = {}
    M = max(e) if e else 0
    for r in range(1, M+1):
        t = e.get(r, ZERO)
        s = ZERO
        for m in range(r, M+1):
            em = e.get(m, ZERO)
            s = s + em.conj_i() * (IQ**(-m)) * Fr(comb(m, r)) * (IQ**(m-r))
        out[r] = t - s
    return out

for p in (2,3):
    for k in (2*p, 2*p+1, 2*p+2):
        ps = Psi_orders(p,p,k)
        high = {r:c for r,c in ps.items() if r>=3 and not c.is_zero()}
        print(f'p={p} k={k}: Psi orders >=3 nonzero: {high}')
    print()
