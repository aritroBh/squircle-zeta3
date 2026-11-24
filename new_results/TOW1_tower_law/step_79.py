from towlib import *
from math import comb
def e_from_data(key, delta, s):
    i,j = key; d = i+j
    vals = [u(i,j,d+t,delta) for t in range(0, s+1)]
    return sum( (vals[t].scale(Fr((-1)**(s-t) * comb(s,t))) for t in range(0, s+1)), ZERO )

for key in [(6,1),(7,0),(5,2),(8,0)]:
    d = sum(key)
    r0 = rho_c(key,0)
    for delta,s in [(6,3),(7,3)]:
        if delta > 2*d-4 or d+s > 14: continue
        e = e_from_data(key, delta, s)
        print(f'{key} d={d}: e({delta},{s})/rho0 = {e/r0}')
    print()
