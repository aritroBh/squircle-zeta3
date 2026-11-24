import importlib, towlib
importlib.reload(towlib)
from towlib import *
bad = []
for k in range(3, 15):
    for d in range(2, min(10,k)+1):
        for key in keys_of_degree(d):
            i,j = key
            if k == d: continue
            for delta in range(0, min(2, 2*d-4)+1):
                lhs = U_e(key,k,delta)
                rhs = Lambda_ij(i,j,k,delta, lambda kk,dp: U_e(kk,k,dp))
                if lhs != rhs:
                    bad.append((k,key,delta,lhs-rhs))
print('T1 for e-built family, all levels 3..14, delta<=2: bad =', len(bad))
for b in bad[:8]: print(b)
