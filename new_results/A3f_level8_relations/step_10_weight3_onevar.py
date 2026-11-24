"""step_10_weight3_onevar.py: proven one-variable Li3 relations at level-8 cross-ratios.

(T3) three-term/Landen, 0<z<1 (theorem: d/dz both sides agree, value at z->0+ is zeta(3)):
  Li3(z) + Li3(1-z) + Li3(z/(z-1)) = zeta(3) + (pi^2/6) log(1-z) - (1/2) log z log^2(1-z) + (1/6) log^3(1-z)
(I3) inversion, z<0 (theorem: d/dz both sides agree, value at z=-1 checks Li3(-1)=-3/4 zeta(3)):
  Li3(z) - Li3(1/z) = -(pi^2/6) log(-z) - (1/6) log^3(-z)
"""
import sys, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from fractions import Fraction as F
from itertools import permutations
from cmzv import POINTS, cross, atom, mul, add
from relations import words, primitive
from qsqrt2 import QS2, A, sunit_log

ZERO, ONE = QS2(0), QS2(1)
LOG2 = {(4,): F(-1)}
LOGA = dict(atom(('log', A)))
def log_atom(z):
    j, m = sunit_log(abs(z))
    return add((LOG2, F(j, 2)), (LOGA, F(m)))
def li3_atom(z): return dict(atom(('li', 3, z)))
def pi2(): return dict(atom(('pi2',)))
def z3(): return dict(atom(('z3',)))

S = {}
for p0, p1, pinf, e in permutations(POINTS, 4):
    if p0 is None or e is None: continue
    if pinf is not None and min(p0, e) < pinf < max(p0, e): continue
    z = cross(e, p0, p1, pinf)
    if z is None or z in (ZERO, ONE): continue
    S.setdefault(z, (p0, p1, pinf, e))

rows, names = [], []
def emit(d, name):
    r = primitive(d)
    if r: rows.append(r); names.append(name)

# T3 instances: need z in (0,1), 1-z in S, w=z/(z-1) in S (w<0), log atoms for z,1-z
nt = 0
for z in list(S):
    if not (ZERO < z < ONE): continue
    if (ONE - z) not in S: continue
    w = z / (z - ONE)
    if w not in S: continue
    try:
        lz = log_atom(z); l1z = log_atom(ONE - z)
    except ValueError:
        continue
    lhs = add((li3_atom(z), 1), (li3_atom(ONE - z), 1), (li3_atom(w), 1))
    rhs = add((z3(), 1), (mul(pi2(), l1z), F(1, 6)),
              (mul(mul(lz, l1z), l1z), F(-1, 2)), (mul(mul(l1z, l1z), l1z), F(1, 6)))
    emit(add((lhs, 1), (rhs, -1)), f'T3 z={z}')
    nt += 1
print(f'T3 instances: {nt}', flush=True)

# I3 instances: z<0, 1/z in S, -z S-unit
ni = 0
for z in list(S):
    if not (z < ZERO): continue
    iz = ONE / z
    if iz not in S: continue
    try:
        lz = log_atom(-z)
    except ValueError:
        continue
    lhs = add((li3_atom(z), 1), (li3_atom(iz), -1))
    rhs = add((mul(pi2(), lz), F(-1, 6)), (mul(mul(lz, lz), lz), F(-1, 6)))
    emit(add((lhs, 1), (rhs, -1)), f'I3 z={z}')
    ni += 1
print(f'I3 instances: {ni}', flush=True)

# numeric check at 60 digits
import check_conversion as CC
from check_conversion import evaluate, approx, mp
worst = mp.mpf(0)
for r, nm in zip(rows, names):
    v = abs(approx(evaluate(dict(r))))
    worst = max(worst, v)
    if v > mp.mpf('1e-40'): print('BAD', nm, v)
print(f'{len(rows)} rows numerically zero, worst={mp.nstr(worst,8)}', flush=True)

# exact rank effect on top of weight3_with_lifts
base = pickle.load(open('weight3_with_lifts.pickle', 'rb'))
cols = base['cols']; known = base['rows']
ids = {w: i for i, w in enumerate(cols)}
from flint import fmpq_mat
def rank(rowlist):
    M = fmpq_mat(len(rowlist), len(cols))
    for i, row in enumerate(rowlist):
        for w, c in row:
            M[i, ids[w]] = int(c)
    _, r = M.rref(inplace=True)
    return r
r0 = rank(known)
print('base rank', r0, 'quotient', len(cols) - r0, flush=True)
for row, nm in zip(rows, names):
    r = rank(known + [row])
    if r > r0:
        print('NEW RANK from', nm)
rall = rank(known + rows)
print(f'with all T3/I3: rank={rall} quotient={len(cols)-rall} (DG 27)', flush=True)
with open('weight3_onevar.pickle', 'wb') as f:
    pickle.dump({'rows': rows, 'names': names}, f)
