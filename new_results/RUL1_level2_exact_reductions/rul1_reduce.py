"""rul1_reduce.py -- exact reduction of word-combos to the named basis, with certificates."""
import pickle
from fractions import Fraction as Fr
import rul1la, rul1_named

_data = None
_rref_cache = {}
def load():
    global _data
    if _data is None:
        with open('rul1_data.pkl','rb') as f: _data = pickle.load(f)
    return _data

def _rref(D, w=None):
    if w is not None and w in _rref_cache:
        return _rref_cache[w]
    rr = [rul1la.Row(dict(d), dict(h), piv) for d, h, piv in D['rref']]
    if w is not None:
        _rref_cache[w] = rr
    return rr

def reduce_combo(w, combo, want_history=False):
    D = load()[w]
    rref = _rref(D, w)
    res, hist = rul1la.reduce_vector(combo, rref)
    free = D['free']
    rv = [res.get(f, Fr(0)) for f in free]
    Minv = D['Minv']
    n = len(free)
    c = [sum(Minv[j][i]*rv[j] for j in range(n)) for i in range(n)]   # c solves rv = M^T c
    fnames = D['fnames']
    coords = {fnames[i]: c[i] for i in range(n) if c[i]}
    if not want_history:
        return coords
    total_h = dict(hist)
    for i, cf in enumerate(c):
        if not cf: continue
        for rid, v in D['hist_named'][i].items():
            total_h[rid] = total_h.get(rid, Fr(0)) - cf*v
    total_h = {k: v for k, v in total_h.items() if v}
    return coords, total_h

def check_certificate(w, combo, coords, hist):
    """exact check: combo - sum coords*named - sum hist*R_id == 0 as word vectors."""
    D = load()[w]
    acc = dict(combo)
    for nm, cf in coords.items():
        combo_nm = dict(rul1_named.NAMED[w])[nm]
        for k, v in combo_nm.items():
            acc[k] = acc.get(k, Fr(0)) - cf*v
    rels = D['rels']
    for rid, v in hist.items():
        fam, tag, r = rels[rid]
        for k, x in r.items():
            acc[k] = acc.get(k, Fr(0)) - v*x
    acc = {k: x for k, x in acc.items() if x}
    return acc == {}
