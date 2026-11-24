import sys, os
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from math import comb

def e_from_data(key, delta, s):
    i,j = key; d = i+j
    base = (2*i+1) if i==j else d
    return sum( u(i,j,base+t, delta).scale(Fr((-1)**(s-t)*comb(s,t))) for t in range(0,s+1) )

def catalan(n):
    return comb(2*n,n)//(n+1)

for s in range(1,7):
    delta = 2*s
    print(f'=== top e({delta},{s})/e00 vs Catalan guess, C_s={catalan(s)}')
    for d in range(s+2, 12):
        line = f'  d={d}: '
        for key in keys_of_degree(d):
            i,j = key
            base = (2*i+1) if i==j else d
            if delta > 2*d-4 or base+s > 14: continue
            e = e_from_data(key, delta, s)
            e00 = e_from_data(key,0,0)
            r = e/e00
            prod = 1
            for jj in range(1,s+1): prod *= (2*d-2*jj-1)
            conj = Fr((-1)**s * catalan(s), prod)
            ok = (r == QSI(conj))
            line += f'{key}:{r}{"" if ok else " CONJ="+str(conj)+" **"}  '
        print(line)
