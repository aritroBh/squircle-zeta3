import importlib, towlib
importlib.reload(towlib)
from towlib import *
bad2 = check_newton(9)
print('Newton identities d<=9, delta<=2: bad =', len(bad2))
for b in bad2[:8]:
    key,delta,sigma,lhs,rhs = b
    print(f'{key} d={delta} s={sigma}: lhs-rhs = {lhs-rhs}')
