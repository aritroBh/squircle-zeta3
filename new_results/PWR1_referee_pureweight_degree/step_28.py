import pwr1_engine as E, sympy as sp
for n in [3,4,5,6]:
    tot, pure, B = E.total_log(n)
    gam_ok = all(sp.expand(tot.coeff(k)).coeff(E.GAM,1)==0 for k in range(E.NMAX+1))
    w_ok = all(E.weight_check(pure.coeff(k), k) for k in range(E.NMAX+1))
    w1 = pure.coeff(1)
    print(f"n={n}: gamma-free: {gam_ok}, pure weight: {w_ok}, w1={w1}")
    print(f"   w2 = {pure.coeff(2)}")
    print(f"   w3 = {pure.coeff(3)}")
