import sys, os
from fractions import Fraction as Fr
from functools import lru_cache
sys.path.insert(0, "<project>/research/10_frontier/W1_corner_family/B_taylor_pi4")
from hpl import coeffs_exact

@lru_cache(maxsize=None)
def my_hpl_coeff(w, m):
    w = tuple(w)
    if m <= 0:
        return Fr(0)
    if len(w) == 1:
        a = w[0]
        assert a in (1, -1)
        return Fr(a**(m - 1), m)
    a, rest = w[0], w[1:]
    if a == 0:
        return my_hpl_coeff(rest, m)/m
    return sum((a**(m - 1 - i))*my_hpl_coeff(rest, i) for i in range(1, m)) / m

fails = []
WORDS = [(1,), (-1,), (0,1), (0,-1), (1,1), (1,-1), (-1,1), (-1,-1), (0,0,1), (1,0,1), (0,1,-1), (1,-1,1), (0,0,0,-1), (1,0,-1,1), (-1,-1,0,1,)]
for w in WORDS:
    ce = coeffs_exact(w, 40)
    for m in range(0, 40):
        if my_hpl_coeff(w, m) != ce[m]:
            fails.append((w, m, my_hpl_coeff(w, m), ce[m]))
print("fails:", fails[:10])
