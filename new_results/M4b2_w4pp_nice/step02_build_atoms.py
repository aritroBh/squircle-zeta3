"""M4b2 step02: classical weight-4 basis atoms B_k as exact level-8 word vectors (NICE atom set).
Atoms: the SAME 22 weight-4 products as M4b (products of {log2, L}, {pi^2, Li2(b)},
{z3, Li3(b), Li3(-b), Li3(1/sqrt2)}), plus Li4 at EXACTLY the seven Theorem-W4 arguments
  {b, -b, -b^2, 1/sqrt2, 1/2, 1-1/sqrt2, 2-sqrt2}   (b = sqrt2-1),
in this order. Each Li4 argument is checked to convert via cmzv.find_cross (cross-ratio of the
8 real points) and to build via cmzv.atom; any failure is reported. Vectors via cmzv.atom/mul.
Saves stage2_atoms_nice.pkl."""
import sys, json, pickle, time, itertools
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/research/08_theory/M4b2_w4pp_nice")
import m4b2_common as mc
from fractions import Fraction as F
from qsqrt2 import QS2 as Q, SQ2 as R, A, B
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
mc.log(f'{len(prods)} weight-4 product atoms (same as M4b)')

# the seven NICE Li4 arguments, task order: b, -b, -b^2, 1/sqrt2, 1/2, 1-1/sqrt2, 2-sqrt2
seven = [B, -B, -B*B, R/2, Q(F(1, 2)), Q(1) - R/2, Q(2) - R]
mc.log(f'seven Li4 arguments: {[str(z) for z in seven]}')

# every one must convert via cmzv.find_cross and build via cmzv.atom
fails = []
for z in seven:
    try:
        cr = cm.find_cross(z)
        mc.log(f'  Li4({z}): find_cross OK, cross=(p0,p1,pinf,e)={tuple(str(x) for x in cr)}')
    except ValueError as e:
        fails.append((str(z), str(e)))
        mc.log(f'  Li4({z}): find_cross FAILED: {e}')
mc.log(f'find_cross converts {len(seven) - len(fails)}/{len(seven)} of the seven NICE arguments')
assert not fails, fails

atoms = []
for name, ats in prods:
    v = {(): F(1)}
    for a in ats:
        v = cm.mul(v, cm.atom(a))
    atoms.append((name, [[a[0]] + [str(x) for x in a[1:]] if a[0] in ('log', 'li') else [a[0]] for a in ats], v))
for z in seven:
    atoms.append((f'Li4({z})', [['li', 4, str(z)]], cm.atom(('li', 4, z))))
cols = set(mc.words4())
for name, _, v in atoms:
    assert all(w in cols for w in v), ('non-admissible word', name)
    assert all(len(w) == 4 for w in v), ('wrong weight', name)
mc.log(f'{len(atoms)} atoms (22 products + 7 Li4), all vectors admissible weight-4. s={time.time()-t0:.1f}')
with open(str(mc.HERE/'stage2_atoms_nice.pkl'), 'wb') as f:
    pickle.dump({'atoms': [(name, enc, v) for name, enc, v in atoms],
                 'seven': [str(z) for z in seven]}, f)
mc.log('saved stage2_atoms_nice.pkl')
