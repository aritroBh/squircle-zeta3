from towlib import *
for d in range(3, 10):
    for key in keys_of_degree(d):
        if key[0]==key[1]: continue
        w = e_off(key,2,1)
        r0 = rho_c(key,0)
        # ratio w / rho(0)?
        print(f'{key}: w = {w}   rho0 = {r0}   w/rho0 = {w/r0 if not r0.is_zero() else None}')
    print()
