import sys, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from itertools import permutations
from cmzv import POINTS, cross
from qsqrt2 import QS2, A, B, SQ2, sunit_log

ZERO, ONE = QS2(0), QS2(1)
# real cross ratios of the 8-point configuration (as in step_10)
S = {}
for p0, p1, pinf, e in permutations(POINTS, 4):
    if p0 is None or e is None: continue
    if pinf is not None and min(p0, e) < pinf < max(p0, e): continue
    z = cross(e, p0, p1, pinf)
    if z is None or z in (ZERO, ONE): continue
    S.setdefault(z, (p0, p1, pinf, e))
print('|S| real cross-ratios =', len(S))
Ssorted = sorted(S, key=float)
for z in Ssorted:
    try:
        j, m = sunit_log(z)
        tag = f'2^({j}/2)*a^{m}'
    except ValueError:
        tag = 'NOT S-unit'
    print(f'{float(z):+.8f}  {z}   {tag}')
