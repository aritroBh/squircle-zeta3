import sys
wd = '<project>/research/08_theory/TOW5R_review'
sys.path.insert(0, wd+'/tow5_copy')
import sympy as sp
import tow2lib as T2, tow5lib as T5
d, m = T2.d, T2.m
i = sp.I
A9 = sp.cancel(T5.arho(9)); A10 = sp.cancel(T5.arho(10))
n9, _ = sp.fraction(A9); n10, _ = sp.fraction(A10)
P9 = sp.Poly(sp.expand(n9), d, m); P10 = sp.Poly(sp.expand(n10), d, m)
t9 = sorted(P9.terms(), key=lambda kv: -kv[0][0])[:4]
t10 = sorted(P10.terms(), key=lambda kv: -kv[0][0])[:4]
print('N9 top terms (monom, coef):')
for mn, c in t9: print('  d^%d m^%d : %s' % (mn[0], mn[1], c))
print('N10 top terms:')
for mn, c in t10: print('  d^%d m^%d : %s' % (mn[0], mn[1], c))
# claimed: N_9 = -256 i d^13 + 12288 i d^12 - 4608 d^11 - 273920 i d^11
# claimed: N_10 = -1024 d^15 + 61440 d^14 + (-1726720 + 23040 i) d^13
