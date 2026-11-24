from fractions import Fraction as Fr
from itertools import permutations
INF = None
def cr(x, p0, p1, pinf):
    # cross(x; p0,p1,pinf) = (x-p0)(p1-pinf)/((x-pinf)(p1-p0)), mapping p0->0,p1->1,pinf->inf
    def sub(a,b):
        if a is INF and b is INF: return Fr(0)
        if a is INF: return INF
        if b is INF: return INF
        return a-b
    def div(a,b):
        if a is INF and b is INF: return Fr(1)
        if a is INF: return INF
        if b is INF: return Fr(0)
        return a/b
    num = sub(x,p0); den = sub(x,pinf)
    f1 = div(num, den)
    f2 = div(sub(p1,pinf), sub(p1,p0))
    if f1 is INF or f2 is INF: return INF
    return f1*f2

PTS = [Fr(0), Fr(1,2), Fr(1), Fr(-1), INF]
names = {Fr(0):'0', Fr(1,2):'1/2', Fr(1):'1', Fr(-1):'-1', INF:'inf'}
vals = {}
for p0,p1,pinf,e in permutations(PTS,4):
    v = cr(e, p0, p1, pinf)
    vals.setdefault(v, []).append((e,p0,p1,pinf))
targets = [Fr(1,3), Fr(-1,3), Fr(2,3), Fr(1,2), Fr(-1,2), Fr(-2), Fr(1,9), Fr(3), Fr(2), Fr(-1), Fr(4,3), Fr(3,2), Fr(1,4)]
for t in targets:
    L = vals.get(t, [])
    print(t, ':', len(L), 'representations;', 'e.g.' , [tuple(names[q] for q in rep) for rep in L[:2]] if L else '')
print()
print('all cross-ratio values:', sorted([v for v in vals if v is not None], key=str))
print('total distinct:', len(vals))
