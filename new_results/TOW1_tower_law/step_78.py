from towlib import *
def e_from_data(key, delta, s):
    """Newton coefficient e(delta,s) from tower data: fit u(d..d+s)."""
    i,j = key; d = i+j
    # solve for e's from values u(d),...,u(d+s)
    from fractions import Fraction as Fr
    vals = [u(i,j,d+t,delta) for t in range(0, s+1)]
    # forward differences in binom basis: e(s) = sum_t (-1)^{s-t} binom(s,t) vals[t]... 
    # since u(d+t) = sum_r binom(t,r) e(r): e(s) = sum_t (-1)^{s-t} binom(s,t) u(d+t)
    from math import comb
    return sum( (vals[t].scale(Fr((-1)**(s-t) * comb(s,t))) for t in range(0, s+1)), ZERO )

for key in [(4,1),(5,0),(5,2),(6,1),(7,0),(8,0),(6,3),(7,2)]:
    d = sum(key)
    if key[0]==key[1]: continue
    r0 = rho_c(key,0)
    out = []
    for delta,s in [(2,1),(3,1),(4,2),(5,2)]:
        if delta > 2*d-4: continue
        try:
            e = e_from_data(key, delta, s)
            out.append(f'e({delta},{s})={e}  ratio={e/r0}')
        except Exception as ex:
            out.append(f'e({delta},{s}): {ex}')
    print(key, out)
