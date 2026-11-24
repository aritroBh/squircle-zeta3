"""r07: claims 5 and C1. For the 38 used frame-w5 rows and for the frame-w5 rows that equal u.f literally (r05),
regenerate from stored provenance with frame_ds.identity/convert and the referee's own domain checks
(every factor argument real with 0<|z|<1; cumulative products and all vertex products are 1/t for frame points |t|>1).
For the u=3,5 conjugate provenance (sigma on Q(sqrt2) applied to the factor arguments) report the largest |z|."""
from rcommon import *
import pickle
from functools import reduce
from itertools import product
L = 'r07.log'
sys.path.insert(0, T + '/A3h4_frame_ds'); sys.path.insert(0, T + '/A3c_claude_elementary/Q')
import frame_ds as fd
from cmzv import POINTS, cross
from qsqrt2 import QS2 as Q
VW5 = pickle.load(open(T + '/W5C1_certificate_scoping/validated_w5.pkl', 'rb'))
X = pickle.load(open(HERE + '/r05_X.pkl', 'rb'))
def qf(x): return float(x.p) + float(x.q) * 2 ** .5
def domain(g, fs):
    p0, p1, pinf = g
    ts = {cross(p, p0, p1, pinf): p for p in POINTS}
    bs = {1 / t for t in ts if t is not None and abs(qf(t)) > 1}
    zs = [z for f in fs for k, z in f]
    if not all(0 < abs(qf(z)) < 1 for z in zs): return None, 'argument with |z|>=1'
    cums = []
    for f in fs:
        b = Q(1); s = []
        for k, z in f: b = b * z; s.append(b)
        if not all(x in bs for x in s): return None, 'cumulative product not 1/t, |t|>1'
        cums.append(s)
    for vertex in product(*((Q(1),) + tuple(s) for s in cums)):
        if any(x != Q(1) for x in vertex) and reduce(lambda a, b: a * b, vertex, Q(1)) not in bs: return None, 'vertex product'
    return ts, ''
used = sorted({i for i, u in X['nonlit']} | {i for i, u in X['lit']})
check = sorted(set(used) | set(X['lit'].values()))
mz = 0; bad = []
for i in check:
    g, fs = VW5['provenance'][i]
    ts, why = domain(g, fs)
    if ts is None: bad.append((i, why)); continue
    raw = fd.primitive(fd.convert(g, ts, fd.identity(fs)))
    if raw != VW5['rows'][i]: bad.append((i, 'regenerated row differs'))
    mz = max(mz, max(abs(qf(z)) for f in fs for k, z in f))
log(f'{len(used)} used frame-w5 rows + literal-image rows {sorted(set(X["lit"].values()))}: regenerated from provenance with domain checks; failures {bad}; max |z| {mz:.4f}', L)
# conjugate arguments (sigma) of the non-literal rows
cz = []
for i in sorted({i for i, u in X['nonlit']}):
    g, fs = VW5['provenance'][i]
    cz.append((i, max(abs(qf(z.conj())) for f in fs for k, z in f)))
log(f'non-literal rows: max |sigma(z)| per row: {[(i, round(v, 3)) for i, v in cz]}; all > 1: {all(v > 1 for _, v in cz)}', L)
