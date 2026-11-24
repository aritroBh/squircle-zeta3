"""W4RX part (b): independent regeneration of the relation rows and comparison with the rows
the certificate uses.  Reads as data: D_certificate.json (labels), weight{1,2}.pickle (the lower
bases: lists of integer rows), validated_w{2,3,4}.pickle (frame provenance), and
data/original_rows.json (output of the rerun of the original code; comparison only)."""
import sys, json, pickle, time
sys.dont_write_bytecode = True
from w4rx_rel import *
T0 = time.time()
def log(s): print(f'[{time.time()-T0:6.1f}s] {s}', flush=True)
ROOT = '<project>/research/08_theory'
orig = json.load(open('data/original_rows.json'))
def orow(lst): return tuple((tuple(w), int(c)) for w, c in lst)
cert = json.load(open(ROOT + '/A3h5_exact_certificate/D_certificate.json'))

class Stub:                      # stands in for qsqrt2.QS2 when unpickling data (no import)
    def __setstate__(self, st):
        st = st[1] if isinstance(st, tuple) else st
        self.k = K(st['p'], st['q'])
class DataUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == 'qsqrt2' and name == 'QS2': return Stub
        if module in ('fractions', 'collections', 'builtins') or module.startswith('copyreg'):
            return super().find_class(module, name)
        raise pickle.UnpicklingError(f'refused class {module}.{name}')
def kk(x): return None if x is None else x.k
def load_data(path): return DataUnpickler(open(path, 'rb')).load()

# ---------------------------------------------------------------- standard rows, weights 1-4
res = {}
lower_sel = {}
for n in (1, 2, 3):
    rows, fam = standard_rows(n, dict(lower_sel))
    cols, oldrows, basis = load_data(ROOT + f'/A3c_claude_elementary/Q/weight{n}.pickle')
    basis = [tuple((tuple(w), int(c)) for w, c in r) for r in basis]
    mine = set(rows)
    inset = all(r in mine for r in basis)
    same_all = [tuple(r) for r in rows] == [tuple((tuple(w), int(c)) for w, c in r) for r in oldrows]
    log(f'weight {n}: regenerated {len(rows)} rows (families: {Counter(f[0] for f in fam)}); '
        f'identical ordered list to the notes\' generated rows: {same_all}; '
        f'all {len(basis)} rows of the notes\' selected basis are regenerated relations: {inset}')
    assert inset
    res[f'w{n}'] = dict(rows=len(rows), same_list=same_all, basis=len(basis), basis_in=inset)
    lower_sel[n] = basis
std4, fam4 = standard_rows(4, {})
same4 = std4 == [orow(r) for r in orig['std_all']]
log(f'weight 4: regenerated {len(std4)} standard rows ({Counter(f[0] for f in fam4)}); '
    f'identical ordered list to the original code\'s: {same4}')
res['w4'] = dict(rows=len(std4), same_list=same4)

# ---------------------------------------------------------------- rows used by the certificate
def frame_row(n, i, extra_word, table):
    data = load_data(ROOT + f'/A3h4_frame_ds/validated_w{n}.pickle')
    g, fs = data['provenance'][i]
    p0, p1, pinf = map(kk, g)
    factors = tuple(tuple((s, z.k) for s, z in f) for f in fs)
    ts = {}
    for p in CONFIG:
        t = cross(p, p0, p1, pinf)
        ts[repr(t)] = (t, p)
    assert len(ts) == 8 and cross(p0, p0, p1, pinf) == K(0) and cross(p1, p0, p1, pinf) == K(1) and cross(pinf, p0, p1, pinf) is None
    letters_ok = {repr(t): p for t, p in ts.values() if t is not None and abs(t) > K(1)}
    # absolute convergence and letters: every cumulative product b of every factor and of
    # every stuffle term is real with |b|<1, and 1/b is the t-value of a configuration point.
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
    # identity: sum_t (-1)^{dt} G(t) [stuffle side, = product of series]
    #         - (-1)^{sum d} G(u1) sh G(u2) sh ...  [shuffle side, = product of integrals]
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
    # conversion: dt/(t-a) = dlog((x-g(a))/(x-pinf)), path x: p0 -> p1 (image of t in [0,1])
    g_of = {repr(K(0)): p0, **letters_ok}
    for w in ident:
        assert all(repr(a) in g_of for a in w)
        assert w[-1] != K(0) and w[0] != K(1)
    # homotopy check in the t coordinate: no singular point of any form has t in (0,1)
    for t, p in ts.values():
        if t is not None and K(0) < t < K(1):
            assert repr(t) not in letters_ok
    out = {}
    for w, c in ident.items():
        forms = [pullback(g_of[repr(a)], pinf) for a in w]
        out = vadd((out, 1), (path_integral(forms, CAY[p0], CAY[p1]), c))
    raw = primitive_H(to_H(out))
    assert raw == tuple((tuple(w), int(c)) for w, c in data['rows'][i]), 'frame row differs from stored row'
    lifted = primitive_H(to_H(vmul(from_H(dict(raw)), {g_of_h(tuple(extra_word))[0]: 1})))
    table.append(dict(n=n, i=i, lift=tuple(extra_word), frame=(p0, p1, pinf), factors=factors,
                      maxmod=maxmod, ts={repr(p): repr(t) for t, p in ts.values()},
                      ident_terms=len(ident)))
    return lifted

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

my_rows = []; table = []; match = Counter(); fams = Counter()
for label, orig_row in zip(cert['labels'], orig['rows']):
    if label[0] == 'T0':
        lab = label[1]
        if lab[0] == 'std':
            row = std4[lab[1]]; fams['std:' + fam4[lab[1]][0]] += 1
        else:
            k = {'lift1': 1, 'lift2old': 2}[lab[0]]
            word = tuple(lab[2])
            assert len(word) == 4 - k and word[0] != 0 and word[-1] != -1, 'lift word not convergent'
            row = primitive_H(to_H(vmul(from_H(dict(lower_sel[k][lab[1]])), {g_of_h(word)[0]: 1})))
            fams[lab[0]] += 1
    else:
        _, n, i, w = label
        assert len(w) == 4 - n and (not w or (w[0] != 0 and w[-1] != -1))
        row = frame_row(n, i, w, table); fams[f'frame_w{n}'] += 1
    my_rows.append(row)
    match[row == orow(orig_row)] += 1
log(f'certificate rows rebuilt independently: {dict(fams)}')
log(f'exact agreement with the rows of the original code: {dict(match)}')
res['cert_rows'] = dict(fams); res['agree'] = dict(match)
pickle.dump({'rows': my_rows, 'labels': cert['labels'], 'frame_table': table, 'res': res,
             'std_fam_used': [fam4[l[1][1]] for l in cert['labels'] if l[0] == 'T0' and l[1][0] == 'std']},
            open('data/check_b.pickle', 'wb'))
log('saved data/check_b.pickle')
