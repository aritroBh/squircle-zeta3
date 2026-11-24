"""Regenerate frame rows from their stored provenance (frame (p0,p1,pinf) and stuffle factors, read as DATA)
with the conversion code of the earlier referee W4RX (w4rx_core: own Q(sqrt2), Cayley map, path conversion;
imports no module of the notes). frame_row() below is W4RX check_b.py's frame_row, copied with three changes:
it takes the (provenance, stored row) pair instead of a file name, it does not lift, and failed assertions
are returned as a reason instead of raised. Covers all 288 validated weight-5 frame rows (W5C1) and all
A3h4 validated rows of weights 2, 3, 4. Every domain/convergence assertion of W4RX is kept."""
from common import *
sys.path.insert(0, T + '/W4RX_independent_check')
from w4rx_core import *
from collections import Counter
from functools import lru_cache
LOG = 'regen_frames.log'

class Stub:
    def __setstate__(self, st):
        st = st[1] if isinstance(st, tuple) else st
        self.k = K(st['p'], st['q'])
class KUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == 'qsqrt2' and name == 'QS2': return Stub
        if module in ('fractions', 'collections', 'builtins') or module.startswith('copyreg'):
            return super().find_class(module, name)
        raise pickle.UnpicklingError(f'refused class {module}.{name}')
def kk(x): return None if x is None else x.k

@lru_cache(maxsize=None)
def kstuffle(u, v):
    if not u: return {v: 1}
    if not v: return {u: 1}
    out = Counter()
    for w, c in kstuffle(u[1:], v).items(): out[(u[0],) + w] += c
    for w, c in kstuffle(u, v[1:]).items(): out[(v[0],) + w] += c
    m = (u[0][0] + v[0][0], u[0][1] * v[0][1])
    for w, c in kstuffle(u[1:], v[1:]).items(): out[(m,) + w] += c
    return dict(out)

def frame_row(prov):
    g, fs = prov
    p0, p1, pinf = map(kk, g)
    factors = tuple(tuple((s, z.k) for s, z in f) for f in fs)
    ts = {}
    for p in CONFIG:
        t = cross(p, p0, p1, pinf)
        ts[repr(t)] = (t, p)
    assert len(ts) == 8 and cross(p0, p0, p1, pinf) == K(0) and cross(p1, p0, p1, pinf) == K(1) and cross(pinf, p0, p1, pinf) is None
    letters_ok = {repr(t): p for t, p in ts.values() if t is not None and abs(t) > K(1)}
    def cum(f):
        out = []; b = K(1)
        for s, z in f: b = b * z; out.append(b)
        return out
    maxmod = K(0)
    for f in factors:
        for s, z in f: assert K(0) < abs(z) < K(1), ('|z|<1 fails', z)
    def st_series(us):
        out = {(): 1}
        for u in us:
            new = Counter()
            for v, c in out.items():
                for w, d in kstuffle(v, u).items(): new[w] += c * d
            out = {w: c for w, c in new.items() if c}
        return out
    prod_st = st_series(factors)
    for t in prod_st:
        for b in cum(t):
            assert abs(b) < K(1) and repr(1 / b) in letters_ok, ('stuffle vertex', b)
            maxmod = max(maxmod, abs(b), key=lambda x: x.mpf())
    def gword(f):
        out = []
        for s, b in zip([s for s, z in f], cum(f)): out.extend([K(0)] * (s - 1) + [1 / b])
        return tuple(out)
    ident = Counter()
    for t, c in prod_st.items(): ident[gword(t)] += c * (-1) ** len(t)
    sf = {(): 1}
    for f in factors:
        new = Counter()
        for v, c in sf.items():
            for w, d in sh(v, gword(f)).items(): new[w] += c * d
        sf = new
    depth = sum(len(f) for f in factors)
    for w, c in sf.items(): ident[w] -= c * (-1) ** depth
    ident = {w: c for w, c in ident.items() if c}
    g_of = {repr(K(0)): p0, **letters_ok}
    for w in ident:
        assert all(repr(a) in g_of for a in w)
        assert w[-1] != K(0) and w[0] != K(1)
    for t, p in ts.values():
        if t is not None and K(0) < t < K(1):
            assert repr(t) not in letters_ok
    out = {}
    for w, c in ident.items():
        forms = [pullback(g_of[repr(a)], pinf) for a in w]
        out = vadd((out, 1), (path_integral(forms, CAY[p0], CAY[p1]), c))
    return primitive_H(to_H(out)), float(maxmod.mpf())

def run(label, path, idxs=None):
    d = KUnpickler(open(path, 'rb')).load()
    rows, prov = d['rows'], d['provenance']
    idxs = range(len(rows)) if idxs is None else idxs
    res = {}
    t0 = time.time()
    for i in idxs:
        stored = tuple((tuple(w), int(c)) for w, c in rows[i])
        try:
            r, mm = frame_row(prov[i])
            res[i] = ('equal' if r == stored else 'DIFFERS', mm)
        except AssertionError as ex:
            res[i] = ('ASSERT ' + repr(ex.args)[:120], None)
    c = Counter(v[0] for v in res.values())
    mm = max((v[1] for v in res.values() if v[1] is not None), default=None)
    log(f'{label}: {len(res)} rows regenerated from provenance; outcome {dict(c)}; max |stuffle vertex| {mm}; ({time.time()-t0:.0f}s)', LOG)
    bad = {i: v for i, v in res.items() if v[0] != 'equal'}
    if bad: log(f'{label}: not equal: {bad}', LOG)
    return res

if __name__ == '__main__':
    out = {}
    for n in (2, 3, 4):
        out[n] = run(f'A3h4 validated_w{n}', T + f'/A3h4_frame_ds/validated_w{n}.pickle')
    out[5] = run('W5C1 validated_w5 (framew5)', T + '/W5C1_certificate_scoping/validated_w5.pkl')
    pickle.dump(out, open(HERE + '/regen_frames.pkl', 'wb'))
