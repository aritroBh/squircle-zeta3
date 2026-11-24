"""r04: every frame row used (66 weight-5 rows in the frame bases, 13 A3h4 rows of weight 2-4 lifted in the frame bases
or in the weight-4 span) regenerated from its stored provenance with A3h4's frame_ds.identity/convert (REUSED project
code: the conversion is hypothesis (b)) and with MY OWN domain check: every factor argument real with 0<|x|<1, every
cumulative product and every nonempty vertex product of cumulative products equal to 1/t for a frame point t with |t|>1.
Also: t06's candidates() test is equivalent (checked on these rows), and the u=3,5 images of the 37 rows are not in domain."""
from rc import *
from functools import reduce
from itertools import product
sys.path.insert(0, T + '/A3h4_frame_ds'); sys.path.insert(0, T + '/A3c_claude_elementary/Q')
import frame_ds as fd
from cmzv import POINTS, cross
from qsqrt2 import QS2 as Q
LOG = 'r04.log'
def qf(x): return float(x.p) + float(x.q) * 2 ** .5
def domain(g, fs):
    p0, p1, pinf = g
    ts = {cross(p, p0, p1, pinf): p for p in POINTS}
    inv = {1 / t for t in ts if t is not None and abs(qf(t)) > 1}
    for f in fs:
        for k, x in f:
            if not (0 < abs(qf(x)) < 1): return None, 'argument'
    cums = []
    for f in fs:
        b = Q(1); s = []
        for k, x in f: b = b * x; s.append(b)
        if not all(y in inv for y in s): return None, 'cumulative'
        if not all(abs(qf(s[i])) < abs(qf(s[i - 1])) for i in range(1, len(s))): return None, 'not decreasing'
        cums.append(s)
    for vertex in product(*((Q(1),) + tuple(s) for s in cums)):
        if any(y != Q(1) for y in vertex) and reduce(lambda a, b: a * b, vertex, Q(1)) not in inv: return None, 'vertex'
    return ts, ''
cs = pickle.load(open(T + '/W5G2_frame_exact/cert_summary.pkl', 'rb'))['summary']
used5 = sorted({lab[1] for r in cs for lab, a in r['frame_terms'] if lab[0] == 'framew5'})
FR13 = [(2, 3), (3, 3), (3, 5), (3, 6), (4, 4), (4, 6), (4, 8), (4, 9), (4, 11), (4, 12), (4, 40), (4, 41), (4, 43)]
w4u = pickle.load(open(T + '/W5G2_frame_exact/w4cert_summary.pkl', 'rb'))['frame_rows_used']
log(f'frame-w5 rows with nonzero a_j in some certificate: {len(used5)}; A3h4 rows used in weight-4 span (w4cert_summary): {w4u}; '
    f'subset of the 13 named rows: {set(map(tuple, w4u)) <= set(FR13)}', LOG)
V5 = pickle.load(open(T + '/W5C1_certificate_scoping/validated_w5.pkl', 'rb'))
bad = []; mz = 0
for i in used5:
    g, fs = V5['provenance'][i]
    ts, why = domain(g, fs)
    if ts is None: bad.append((5, i, why)); continue
    if fd.primitive(fd.convert(g, ts, fd.identity(fs))) != V5['rows'][i]: bad.append((5, i, 'differs'))
    mz = max(mz, max(abs(qf(x)) for f in fs for k, x in f))
for n, i in FR13:
    d = pickle.load(open(T + f'/A3h4_frame_ds/validated_w{n}.pickle', 'rb'))
    g, fs = d['provenance'][i]
    ts, why = domain(g, fs)
    if ts is None: bad.append((n, i, why)); continue
    if fd.primitive(fd.convert(g, ts, fd.identity(fs))) != d['rows'][i]: bad.append((n, i, 'differs'))
    mz = max(mz, max(abs(qf(x)) for f in fs for k, x in f))
log(f'regenerated {len(used5)} weight-5 + {len(FR13)} weight-2..4 frame rows with own domain check: failures {bad}; max |x| {mz:.4f}', LOG)
# conjugate provenance of the 37 open rows: out of domain?
OPEN = [4, 5, 7, 9, 10, 12, 14, 15, 16, 17, 52, 53, 54, 56, 58, 59, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 87, 91, 95, 103, 104, 107, 113]
cz = [max(abs(qf(x.conj())) for f in V5['provenance'][i][1] for k, x in f) for i in OPEN]
log(f'37 open rows: min over rows of max |sigma(x)| = {min(cz):.3f} (>1: conjugate provenance is outside the domain)', LOG)
