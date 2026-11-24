import sys, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from itertools import permutations
from cmzv import POINTS, cross
from qsqrt2 import QS2, A, B, SQ2, sunit_log

ZERO, ONE = QS2(0), QS2(1)
S = {}
for p0, p1, pinf, e in permutations(POINTS, 4):
    if p0 is None or e is None: continue
    if pinf is not None and min(p0, e) < pinf < max(p0, e): continue
    z = cross(e, p0, p1, pinf)
    if z is None or z in (ZERO, ONE): continue
    S.setdefault(z, (p0, p1, pinf, e))

# full 20-term Kummer argument list with exponents (a,b,c,d) for x^a y^b xi^c eta^d and K4 coef d_i
K20 = [
    (-1,2,1,-1,1, 1), (-1,1,2,1,-1, 1), (1,2,1,-1,-2, 1), (1,1,2,2,-1, 1),
    (1,1,1,0,0, -6), (1,1,1,-1,-1, -6), (-1,1,1,-1,0, -6), (-1,1,1,0,-1, -6),
    (1,1,0,0,1, -3), (1,0,1,1,0, -3), (-1,1,0,-1,1, -3), (-1,0,1,1,-1, -3),
    (-1,1,0,-1,-1, -3), (-1,0,1,-1,-1, -3), (1,1,0,0,-1, -3), (1,0,1,-1,0, -3),
    (1,1,0,0,0, 6), (1,0,1,0,0, 6), (-1,1,0,-1,0, 6), (-1,0,1,0,-1, 6),
]
def mon(s,a,b,c,d,x,y):
    return QS2(s)*x**a*y**b*(ONE-x)**c*(ONE-y)**d

pos = sorted([z for z in S if ZERO < z < ONE], key=float)
print('positive cross-ratios in (0,1):', [str(z) for z in pos])
direct = []
for x in pos:
    for y in pos:
        if x + y > ONE: continue
        args = [mon(s,a,b,c,d,x,y) for s,a,b,c,d,di in K20]
        if all((z in S or z in (ZERO, ONE)) for z in args):
            direct.append((x,y))
print('direct instances (all 20 args in config):', len(direct))
for x,y in direct: print('  ', x, y)
# near miss count
nm = 0
for x in pos:
    for y in pos:
        if x + y > ONE: continue
        args = [mon(s,a,b,c,d,x,y) for s,a,b,c,d,di in K20]
        bad = [z for z in args if z not in S and z not in (ZERO,ONE)]
        if len(bad) <= 4:
            nm += 1
            print(f'near-miss x={x} y={y}: {len(bad)} outside, e.g. {bad[:3]}')
print('near misses (<=4 outside):', nm)
