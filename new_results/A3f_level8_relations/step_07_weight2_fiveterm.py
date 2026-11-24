"""step_07_weight2_fiveterm.py: Rogers five-term instances on real level-8 cross-ratios.

L(z) = Li2(z) + (1/2) log z log(1-z)   (Rogers dilogarithm, real, 0<z<1)
Five-term (Rogers 1907, classical; e.g. Zagier, "The dilogarithm function", eq. (I.1)):
  L(x)+L(y) = L(xy) + L(u) + L(v),  u = x(1-y)/(1-xy), v = y(1-x)/(1-xy),   0<x,y<1.
Every argument that occurs must be a level-8 cross-ratio (atom via cmzv.find_cross)
with z and 1-z S-units in Q(sqrt2) (so log z, log(1-z) are weight-1 MPV combos).
Also reflection: L(z)+L(1-z) = zeta(2) = pi^2/6.
"""
import sys, pickle, time, json
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from fractions import Fraction as F
from itertools import permutations
from cmzv import POINTS, cross, atom, mul, add, regularize
from relations import words, primitive
from qsqrt2 import QS2, A, SQ2, sunit_log

ZERO, ONE = QS2(0), QS2(1)
LOG2 = {(4,): F(-1)}                      # atom(('log',Q(2)))
LOGA = dict(atom(('log', A)))             # log(1+sqrt2)

def log_atom(z):
    """word dict for log z, z>0 an S-unit."""
    j, m = sunit_log(z)
    return add((LOG2, F(j, 2)), (LOGA, F(m)))

def li2_atom(z):
    return dict(atom(('li', 2, z)))

def L(z):
    """Rogers L(z) as a formal weight-2 word dict."""
    return add((li2_atom(z), 1), (mul(log_atom(z), log_atom(ONE - z)), F(1, 2)))

def pi2_atom():
    return dict(atom(('pi2',)))

# --- enumerate admissible cross-ratios S with a valid quadruple
S = {}
for p0, p1, pinf, e in permutations(POINTS, 4):
    if p0 is None or e is None: continue
    if pinf is not None and min(p0, e) < pinf < max(p0, e): continue
    z = cross(e, p0, p1, pinf)
    if z is None or z == ZERO or z == ONE: continue
    if z not in S: S[z] = (p0, p1, pinf, e)
print(f'admissible cross-ratios |S|={len(S)}', flush=True)

GOOD = []
for z in S:
    if ZERO < z < ONE and (ONE - z) in S:
        try:
            sunit_log(z); sunit_log(ONE - z)
        except ValueError:
            continue
        GOOD.append(z)
print(f'GOOD (0,1) cross-ratios with 1-z in S and S-unit logs: {len(GOOD)}')
for z in sorted(GOOD): print('   ', z, '1-z =', ONE - z)

rows = []
names = []
def emit(d, name):
    r = primitive(d)
    if not r:
        print('ZERO ROW:', name); return
    rows.append(r); names.append(name)

# Rogers five-term
n5 = 0
for x in GOOD:
    for y in GOOD:
        xy = x * y
        den = ONE - xy
        u = x * (ONE - y) / den
        v = y * (ONE - x) / den
        for t in (xy, u, v):
            if not (ZERO < t < ONE) or t not in S or (ONE - t) not in S:
                break
        else:
            try:
                for t in (xy, u, v):
                    sunit_log(t); sunit_log(ONE - t)
            except ValueError:
                continue
            d = add((L(x), 1), (L(y), 1), (L(xy), -1), (L(u), -1), (L(v), -1))
            emit(d, f'Rogers5T x={x} y={y}')
            n5 += 1
print(f'Rogers five-term instances: {n5}', flush=True)
# reflection
nr = 0
for z in GOOD:
    if (ONE - z) in GOOD and z <= ONE - z:
        d = add((L(z), 1), (L(ONE - z), 1), (pi2_atom(), F(-1, 6)))
        emit(d, f'reflection z={z}')
        nr += 1
print(f'reflection instances: {nr}', flush=True)

with open('weight2_fiveterm.pickle', 'wb') as f:
    pickle.dump({'rows': rows, 'names': names}, f)

# numeric verification at ~60 digits
import check_conversion as CC
from check_conversion import evaluate, approx, mp
worst = mp.mpf(0)
for r, nm in zip(rows, names):
    v = abs(approx(evaluate(dict(r))))
    worst = max(worst, v)
    if v > mp.mpf('1e-40'): print('BAD', nm, v)
print(f'all {len(rows)} candidate rows numerically zero, worst={mp.nstr(worst,8)}', flush=True)

# exact rank effect
from flint import fmpq_mat
known = pickle.load(open('audit_weight2.pickle', 'rb'))['rows']
cols = words(2)
ids = {w: i for i, w in enumerate(cols)}
def rank_of(rowlist):
    M = fmpq_mat(len(rowlist), len(cols))
    for i, row in enumerate(rowlist):
        for w, c in (row if isinstance(row, tuple) else tuple(row.items())):
            M[i, ids[w]] = int(c) if isinstance(row, tuple) else F(c)
    rr, r = M.rref(inplace=True)
    return rr, r
rr0, r0 = rank_of(known)
allrows = known + rows
rr1, r1 = rank_of(allrows)
print(f'rank known={r0} (quotient {len(cols)-r0}); with five-term rows rank={r1} (quotient {len(cols)-r1})', flush=True)
