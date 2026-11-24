from common import *
import itertools
t0 = time.time()
d89 = pickle.load(open(T+'I_reductions_w89.pkl', 'rb'))
d9 = pickle.load(open(T+'I_reductions_w9.pkl', 'rb'))
basis = {8: d89['basis8'], 9: d89['basis9']}
reds = {8: d89['reds'], 9: d9}

def li_half(n):  # Li_n(1/2) as word combination (the notes' definition; checked numerically in run_num.py)
    out = {}
    for c in itertools.product((0, -1), repeat=n-1):
        w = c + (-1,); out[w] = out.get(w, F(0)) + (-1)**sum(1 for a in c if a == -1)
    return {k: v for k, v in out.items() if v}
PRIM = {'log2': {(-1,): F(1)},
        'X6': zbar(((5, -1), (1, -1))), 'X7a': zbar(((5, -1), (1, 1), (1, 1))), 'X7b': zbar(((4, -1), (2, 1), (1, 1))),
        'X8a': zbar(((6, -1), (1, 1), (1, 1))), 'X8b': zbar(((5, -1), (2, 1), (1, 1))),
        'X8c': zbar(((5, -1), (1, 1), (1, 1), (1, 1))), 'X8d': zbar(((4, -1), (2, 1), (1, 1), (1, 1))),
        'zbar(7,1,1)': zbar(((7, -1), (1, 1), (1, 1))), 'zbar(6,2,1)': zbar(((6, -1), (2, 1), (1, 1))),
        'zbar(5,3,1)': zbar(((5, -1), (3, 1), (1, 1))), 'zbar(5,1,1,1,1)': zbar(((5, -1),)+((1, 1),)*4),
        'zbar(4,2,1,1,1)': zbar(((4, -1), (2, 1), (1, 1), (1, 1), (1, 1))),
        'X9fp': {(-1, -1, -1, 0, 0, 0, -1, 0, -1): F(1)}}
def prim(nm):
    if nm in PRIM: return PRIM[nm]
    if nm[0] == 'z': return zeta_word(int(nm[1:]))
    if nm.startswith('Li'): return li_half(int(nm[2:]))
    raise KeyError(nm)
def name_vec(name): return mono(*[prim(x) for x in name.split('*')])

res = {}
for W in (8, 9):
    data, B = load_rows_meta(W)
    bvec = [name_vec(nm) for nm in basis[W]]
    wt_ok = all(all(len(w) == W for w in v) for v in bvec)
    print(f"weight {W}: basis size {len(bvec)}; all basis words have weight {W}: {wt_ok}", flush=True)
    for m in range(W-1):
        j = W-1-m
        h = norm_h(pickle.load(open(T+f'cert_I_{m}_{j}.pkl', 'rb')), len(B))
        coeffs = reds[W][(m, j)]
        assert len(coeffs) == len(bvec)
        V = {(0,)*(m+1)+(-1,)*j: F(fact(m)*fact(j))}      # I_{m,j} = j! m! G(0^{m+1}(-1)^j)
        for c, bv in zip(coeffs, bvec):
            if c: V = vsub(V, bv, F(c))
        ok, mism = verify_certificate(h, V, data, B, crosscheck=True)
        res[(m, j)] = ok
        print(f"  I_{{{m},{j}}} = sum c_b b   [w={W}] certificate EXACT: {ok} (relation-row mismatches {mism})  t={time.time()-t0:.0f}s", flush=True)
print("all I-reduction certificates (w8,w9) exact:", all(res.values()), len(res))
