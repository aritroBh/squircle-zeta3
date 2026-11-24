"""step_02: classical weight-4 basis atoms B_k as exact level-8 word vectors.
  (a) the 22 weight-4 products of {log2, L}, {pi^2, Li2(b)}, {z3, Li3(b), Li3(-b), Li3(1/sqrt2)}
      as in paper/squircle_note/scripts/m3_pslq.py;
  (b) Li4(z) for every real S-unit z < 1 of Q(sqrt2) with 1-z an S-unit (m3_pslq bounds
      j in [-6,6], m in [-4,4], sign +-), split by whether cmzv.find_cross converts it
      (cross-ratio of the 8 real points, convergence condition on pinf).
Vectors via cmzv.atom and cmzv.mul. Saves stage2_atoms.pkl."""
import sys, json, pickle, time, itertools, functools
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/research/08_theory/M4b_w4pp_exact")
import m4b_common as mc
from fractions import Fraction as F
from qsqrt2 import QS2 as Q, SQ2 as R, A, B, sunit_log
SQ2 = R
import cmzv as cm

t0 = time.time()
log2 = ('log', Q(2)); L = ('log', A)
pi2 = ('pi2',); Li2b = ('li', 2, B)
z3 = ('z3',)
Li3b = ('li', 3, B); Li3mb = ('li', 3, -B); Li3h = ('li', 3, R/2)
w1 = [('log2', log2), ('L', L)]
w2 = [('pi^2', pi2), ('Li2(b)', Li2b)]
w3 = [('z3', z3), ('Li3(b)', Li3b), ('Li3(-b)', Li3mb), ('Li3(1/sqrt2)', Li3h)]

prods = []
for c in itertools.combinations_with_replacement(range(2), 4):
    prods.append(('*'.join(w1[i][0] for i in c), tuple(w1[i][1] for i in c)))
for (n2, v2) in w2:
    for c in itertools.combinations_with_replacement(range(2), 2):
        prods.append((n2 + '*' + '*'.join(w1[i][0] for i in c), (v2,) + tuple(w1[i][1] for i in c)))
for c in itertools.combinations_with_replacement(range(2), 2):
    prods.append(('*'.join(w2[i][0] for i in c), tuple(w2[i][1] for i in c)))
for (n3, v3) in w3:
    for (n1, v1) in w1:
        prods.append((n3 + '*' + n1, (v3, v1)))
assert len(prods) == 22
mc.log(f'{len(prods)} weight-4 product atoms')

# candidate Li4 arguments (m3_pslq bounds)
cands = set()
for j in range(-6, 7):
    for m in range(-4, 5):
        for sg in (1, -1):
            z = Q(sg) * SQ2 ** j * A ** m
            if not (z < Q(1)) or z == Q(0):
                continue
            try:
                sunit_log(Q(1) - z)
            except ValueError:
                continue
            if abs(z) > Q(1):
                continue
            cands.add(z)
def cmpz(x, y):
    ax, ay = abs(x), abs(y)
    kx = 1 if ax > Q(9, 10) else 0
    ky = 1 if ay > Q(9, 10) else 0
    if kx != ky: return kx - ky
    if ax != ay: return -1 if ax > ay else 1     # larger |z| first
    return -1 if x < y else (x > y)
cands = sorted(cands, key=functools.cmp_to_key(cmpz))
mc.log(f'{len(cands)} candidate Li4 arguments: {[str(z) for z in cands]}')

converts, fails = [], []
for z in cands:
    try:
        cr = cm.find_cross(z)
        converts.append((z, cr))
    except ValueError:
        fails.append(z)
mc.log(f'find_cross converts {len(converts)}; does NOT convert {len(fails)}: {[str(z) for z in fails]}')
for z, cr in converts:
    mc.log(f'  Li4({z}): cross=(p0,p1,pinf,e)={tuple(str(x) for x in cr)}')
seven = [B, -B, -B*B, R/2, Q(F(1,2)), Q(1)-R/2, Q(2)-R]
missing = [z for z in seven if z not in {z for z, _ in converts}]
mc.log(f'the 7 W4 closed-form arguments present among converting: {not missing} (missing: {[str(z) for z in missing]})')
assert not missing

# build vectors
atoms = []
for name, ats in prods:
    v = {(): F(1)}
    for a in ats:
        v = cm.mul(v, cm.atom(a))
    atoms.append((name, [list(a) if a[0] != 'li' else [a[0], a[1], str(a[2])] for a in ats], v))
for z, cr in converts:
    atoms.append((f'Li4({z})', [['li', 4, str(z)]], cm.atom(('li', 4, z))))
cols = set(mc.words4())
for name, _, v in atoms:
    assert all(w in cols for w in v), ('non-admissible word', name)
    assert all(len(w) == 4 for w in v), ('wrong weight', name)
mc.log(f'{len(atoms)} atoms, all vectors admissible weight-4; total support rat-side n/a. s={time.time()-t0:.1f}')
with open(str(mc.HERE/'stage2_atoms.pkl'), 'wb') as f:
    pickle.dump({'atoms': [(name, enc, v) for name, enc, v in atoms],
                 'converts': [str(z) for z, _ in converts], 'fails': [str(z) for z in fails]}, f)
mc.log('saved stage2_atoms.pkl')
