import sys, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from itertools import permutations
from fractions import Fraction as F
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

# candidate S-units x in (0,1) with 1-x an S-unit, |j|,|m|<=6 (B_kummer: 22, provably complete)
cand = []
for j in range(-6, 7):
    for m in range(-6, 7):
        for sg in (1, -1):
            z = QS2(sg)*SQ2**j*A**m
            if ZERO < z < ONE:
                try: sunit_log(ONE - z)
                except ValueError: continue
                cand.append(z)
cand = sorted(set(cand), key=float)
print('candidates:', len(cand))

K20 = [
    (-1,2,1,-1,1, 1), (-1,1,2,1,-1, 1), (1,2,1,-1,-2, 1), (1,1,2,2,-1, 1),
    (1,1,1,0,0, -6), (1,1,1,-1,-1, -6), (-1,1,1,-1,0, -6), (-1,1,1,0,-1, -6),
    (1,1,0,0,1, -3), (1,0,1,1,0, -3), (-1,1,0,-1,1, -3), (-1,0,1,1,-1, -3),
    (-1,1,0,-1,-1, -3), (-1,0,1,-1,-1, -3), (1,1,0,0,-1, -3), (1,0,1,-1,0, -3),
    (1,1,0,0,0, 6), (1,0,1,0,0, 6), (-1,1,0,-1,0, 6), (-1,0,1,0,-1, 6),
]
def mon(s,a,b,c,d,x,y):
    return QS2(s)*x**a*y**b*(ONE-x)**c*(ONE-y)**d

def logvec(z):
    """log|z| as {(e2,eA): coef} rational, z an S-unit."""
    j, m = sunit_log(z)
    return {(1,0): F(j,2), (0,1): F(m)}

def pmul(P, Q):
    out = {}
    for (a,b),c in P.items():
        for (d,e),f in Q.items():
            out[(a+d,b+e)] = out.get((a+d,b+e), F(0)) + c*f
    return out

def ppow(P, n):
    out = {(0,0): F(1)}
    for _ in range(n): out = pmul(out, P)
    return out

def l2la(lg):
    """log = (j/2) log2 + m logA -> (F(j,2), m) coefficients."""
    return lg[(1,0)], lg[(0,1)]

# relation: {'li3': {z: QS2coef}, 'rest': {(e2,eA,epi,ez3): QS2coef}} meaning sum = 0
def mk_Rx(x, y):
    xi, eta = ONE-x, ONE-y
    li3, rest = {}, {}
    for s,a,b,c,d,di in K20:
        m = mon(s,a,b,c,d,x,y)
        cf = QS2(di)*(QS2(a)*xi - QS2(c)*x)
        if cf == ZERO: continue
        if m == ONE: rest[(0,0,0,1)] = rest.get((0,0,0,1), ZERO) + cf
        else: li3[m] = li3.get(m, ZERO) + cf
    # + 3x log(xi) log^2(eta)
    lxi, leta = logvec(xi), logvec(eta)
    poly = pmul(lxi, ppow(leta, 2))
    for (e2,eA),c in poly.items():
        rest[(e2,eA,0,0)] = rest.get((e2,eA,0,0), ZERO) + QS2(3)*x*QS2(c)
    return {'name': f'Rx x={x} y={y}', 'li3': li3, 'rest': rest}

def mk_Ry(x, y):
    xi, eta = ONE-x, ONE-y
    li3, rest = {}, {}
    for s,a,b,c,d,di in K20:
        m = mon(s,a,b,c,d,x,y)
        cf = QS2(di)*(QS2(b)*eta - QS2(d)*y)
        if cf == ZERO: continue
        if m == ONE: rest[(0,0,0,1)] = rest.get((0,0,0,1), ZERO) + cf
        else: li3[m] = li3.get(m, ZERO) + cf
    lxi, leta = logvec(xi), logvec(eta)
    poly = pmul(ppow(lxi,2), leta)
    for (e2,eA),c in poly.items():
        rest[(e2,eA,0,0)] = rest.get((e2,eA,0,0), ZERO) + QS2(3)*y*QS2(c)
    return {'name': f'Ry x={x} y={y}', 'li3': li3, 'rest': rest}

pairs = [(x,y) for i,x in enumerate(cand) for y in cand[i:] if x+y <= ONE]
print('instance pairs:', len(pairs))
rels = [mk_Rx(x,y) for x,y in pairs] + [mk_Ry(x,y) for x,y in pairs]
atoms = set()
for r in rels: atoms |= set(r['li3'])
print('distinct Li3 atoms:', len(atoms), '; aux (not in cross-ratio S):', len([z for z in atoms if z not in S]))
# numeric spot check of Rx at one pair (60-digit polylog)
import mpmath as mp
mp.mp.dps = 60
def numval(r):
    tot = mp.mpf(0)
    for z,c in r['li3'].items():
        tot += float(c)*mp.polylog(3, float(z))
    for (e2,eA,epi,ez),c in r['rest'].items():
        tot += float(c)*(mp.log(2)**e2)*(mp.log(1+mp.sqrt(2))**eA)*(mp.pi**2)**epi*(mp.zeta(3)**ez)
    return tot
worst = max(abs(numval(r)) for r in rels[:20])
print('spot numeric check first 20 rows, worst residual:', mp.nstr(worst, 6))
pickle.dump({'cand':cand,'pairs':pairs,'rels':rels}, open('k3_deriv_raw.pickle','wb'))
print('saved k3_deriv_raw.pickle')
