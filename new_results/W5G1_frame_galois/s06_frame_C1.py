"""s06 (questions A and C1).
A : regenerate every used frame-w5 row from its stored provenance (frame, factors) with frame_ds's own
    identity/convert, re-asserting the absolute-convergence domain (every argument real, 0<|z|<1, every
    cumulative/vertex product in bs), and check it equals the stored row.
C1: for u = 5 (sqrt2 -> -sqrt2 on the real points; on Cayley exponents x -> 5x) and u = 3 (= 5 o 7):
    conjugate the provenance; check whether the conjugate frame is one of the 336 frames and whether the
    conjugate factors satisfy the domain conditions there; if so regenerate and compare with u.row."""
from common import *
import frame_ds as fd
from cmzv import POINTS, cross
from qsqrt2 import QS2 as Q
from functools import reduce
VW5 = pickle.load(open(W5C1 + '/validated_w5.pkl', 'rb'))
split = pickle.load(open(HERE + '/split.pkl', 'rb'))
used_idx = sorted({i for (i, u) in split['frame_pairs']})
FR = [{tuple(w): c for w, c in row} for row in VW5['rows']]
def cj(x): return None if x is None else x.conj()
def neg(x): return None if x is None else -x
ALLF = {g: (ts, bs) for g, ts, bs in fd.frames()}
def in_domain(fs, bs):
    """Exactly the conditions enforced by frame_ds.candidates, for given factors (tuples of (k, z))."""
    cum = []
    for f in fs:
        b = Q(1); seq = []
        for k, z in f:
            b = b * z; seq.append(b)
        if not all(x in bs for x in seq): return False, 'cumulative product not a frame point with |t|>1'
        if not all(abs(float(seq[i].p) + float(seq[i].q) * 2 ** .5) < abs(float(seq[i - 1].p) + float(seq[i - 1].q) * 2 ** .5) for i in range(1, len(seq))):
            return False, 'not decreasing'
        cum.append(seq)
    import itertools
    for vertex in itertools.product(*((Q(1),) + tuple(s) for s in cum)):
        if any(x != Q(1) for x in vertex) and reduce(lambda a, b: a * b, vertex, Q(1)) not in bs:
            return False, 'vertex product not in bs'
    return True, ''
def regen(g, fs):
    ts, bs = ALLF[g]
    ok, why = in_domain(fs, bs)
    if not ok: return None, why
    rel = fd.identity(fs)
    assert all(a in ts for w in rel for a in w)
    return primitive(fd.convert(g, ts, rel)), ''
def key(d):
    k, m = None, None
    return primitive({w: F(c) for w, c in d.items()})
resA = []; resC = []
for i in used_idx:
    g, fs = VW5['provenance'][i]
    row = tuple(sorted((tuple(w), c) for w, c in VW5['rows'][i]))
    r, why = regen(g, fs)
    okA = (r == primitive(dict(row)))
    # also confirm every argument real with |z|<1
    zs = [float(z.p) + float(z.q) * 2 ** .5 for f in fs for k, z in f]
    resA.append((i, okA, max(abs(z) for z in zs)))
    for u in (5, 3):
        if u == 5: gc = tuple(cj(x) for x in g); fsc = tuple(tuple((k, z.conj()) for k, z in f) for f in fs)
        else:      gc = tuple(neg(cj(x)) for x in g); fsc = tuple(tuple((k, -z.conj()) for k, z in f) for f in fs)
        # u=3 = (x->-x)o sigma : the frame map t is invariant under x->-x applied to all of p0,p1,pinf, x;
        # factor arguments are cross-ratio data, so they transform by sigma only.
        if u == 3: fsc = tuple(tuple((k, z.conj()) for k, z in f) for f in fs)
        target = primitive(actv(dict(row), u))
        if gc not in ALLF:
            resC.append((i, u, 'conjugate frame not a frame', None)); continue
        r2, why = regen(gc, fsc)
        if r2 is None:
            resC.append((i, u, 'DOMAIN FAILS: ' + why, [[float(z.p) + float(z.q) * 2 ** .5 for k, z in f] for f in fsc])); continue
        resC.append((i, u, 'regenerated; equals u.row: ' + str(r2 == target), None))
log(f'A: regenerated {len(resA)} used frame-w5 rows from provenance with domain checks; all equal stored rows: '
    f'{all(x[1] for x in resA)}; max |argument| = {max(x[2] for x in resA):.4f}', 's06.log')
import collections
log('C1 outcomes: ' + str(collections.Counter((u, s.split(':')[0]) for i, u, s, _ in resC)), 's06.log')
for i, u, s, z in resC:
    log(f'  row {i} u={u}: {s}  conj args={z}', 's06.log')
pickle.dump(dict(A=resA, C1=resC), open(HERE + '/s06_C1.pkl', 'wb'))
