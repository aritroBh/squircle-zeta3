"""Claude's independent check of KSPAN: (1) fast reducer == slow sympy reference on ALL keys at k = 4, 5;
(2) for each log degree d = 1..5 the degree-d kernel vectors K^(k) satisfy Delta^d_k K = 0 exactly for consecutive k
(i.e. polynomial in k of degree <= d-1) and Delta^(d-1) != 0; (3) Q-rank of {K^(k)} = d (as Q-vectors via a+b*sqrt2).
"""
import sys, pickle, time, os
sys.dont_write_bytecode = True
W2 = '<project>/research/08_theory/WT2_weight2_all_k'
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
sys.path.insert(0, W2)
from fractions import Fraction as Fr
import fastred
from fastred import reduce_k_fast, pf_from_sympy
import sympy as sp

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14


def flat(pf):
    """PF -> {pole_repr: (p, q)} for the simple part, asserting it is real and has nothing at +-i."""
    simple, hp, hm, const = pf
    assert not hp and not hm and const.is_zero()
    out = {}
    for pole, c in simple.items():
        assert c.p2 == 0 and c.p3 == 0
        if c.p0 or c.p1:
            out[repr(pole)] = (c.p0, c.p1)
    return out


def load(k):
    return pickle.load(open(os.path.join(W2, f'binteg_{k}.pkl'), 'rb'))


# (1) fast vs slow on all keys
from ibp_reduce import reduce_k
for k in (4, 5):
    t0 = time.time()
    fast = reduce_k_fast(k, integrand_pf=load(k))
    slow = reduce_k(k, verbose=False)
    keys = set(fast) | set(slow)
    bad = []
    for key in keys:
        if key == (0, 0):
            continue
        f = flat(fast[key]) if key in fast else {}
        s = flat(pf_from_sympy(slow[key])) if key in slow else {}
        if f != s:
            bad.append(key)
    print(f'(1) k={k}: fast vs slow on {len(keys)} keys: {"ALL EQUAL" if not bad else "MISMATCH " + str(bad)} [{time.time()-t0:.0f}s]', flush=True)

# (2), (3)
kern = {}
for k in range(2, KMAX + 1):
    kern[k] = reduce_k_fast(k, integrand_pf=load(k))


def vec(k, d):
    out = {}
    for key, pf in kern[k].items():
        if sum(key) == d:
            for pole, pq in flat(pf).items():
                out[(key, pole)] = pq
    return out


def sub(u, w, c=1):
    out = dict(u)
    for x, (p, q) in w.items():
        a, b = out.get(x, (Fr(0), Fr(0)))
        out[x] = (a - c * p, b - c * q)
    return {x: v for x, v in out.items() if v != (0, 0)}


for d in range(1, 6):
    ks = list(range(d + 1, KMAX + 1)) if d > 1 else list(range(2, KMAX + 1))
    V = [vec(k, d) for k in ks]
    diffs = V
    order = 0
    while any(diffs) and len(diffs) > 1:
        diffs = [sub(diffs[i + 1], diffs[i]) for i in range(len(diffs) - 1)]
        order += 1
    # order = first difference order that vanished identically on all available consecutive k
    cols = sorted({x for v in V for x in v}, key=str)
    M = sp.Matrix([[c for x in cols for c in v.get(x, (0, 0))] for v in V])
    print(f'(2,3) d={d}: k={ks[0]}..{ks[-1]}: Delta^{order} K == 0 identically (remaining diffs {len(diffs)}), '
          f'so polynomial degree = {order-1}; Q-rank of the {len(V)} vectors = {M.rank()}', flush=True)
