import importlib, towlib
importlib.reload(towlib)
from towlib import *
p,k = 2,4
lhs = sp.expand(pf_to_sympy(diag_feed_reduction(p,k)))
rhs = diag_feed_formula(p,k)
diff = sp.expand(lhs-rhs)
print('diff =', diff)
# also print each side's structure
print('lhs =', lhs)
