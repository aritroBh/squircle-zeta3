"""check_full.py -- referee re-derivation of the '38939 cells, 0 bad' statement: every data cell
(key, delta<=min(6,2d-4), k in [k0, min(d+8,52)]), k0=d (off-diagonal) / d+1 (diagonal), compared with the
Newton-form prediction from the claimed e-family (off-diagonal) and ehat-family (diagonal).
Also: diagonal initial level k=2p vs rho closed form a_delta(2p,p)."""
import sys, pickle, time
from math import comb
sys.path.insert(0, '.')
import sympy as sp
import check_data as CD
import tow3lib as T3
import prove_E3b56 as PE
Q = CD.Q
t0 = time.time()
H = pickle.load(open(CD.PKL, 'rb')); H.pop('_meta')
d_, m_ = T3.d, T3.m
p_ = PE.p
evE = {}; evH = {}
for dl in range(0, 7):
    for s in range(0, dl//2+1):
        evE[(dl, s)] = CD.make_eval(T3.e(dl, s), (d_, m_))
        evH[(dl, s)] = CD.make_eval(sp.cancel(sp.together(PE.ehat(dl, s))), (p_,))
print('evaluators built', round(time.time()-t0, 1), flush=True)
cells = 0; bad = 0; diag_cells = 0; init_cells = 0; init_bad = 0
bad_ex = []
for (p, q), v in sorted(H.items()):
    d = p+q
    diag = (p == q)
    for dl in range(0, min(6, 2*d-4)+1):
        rho0 = Q.of(H[(p, q)][0][d])
        P = dl//2
        if not diag:
            coeffs = [evE[(dl, s)](d, p) for s in range(P+1)]
            k0 = d
        else:
            coeffs = [evH[(dl, s)](p) for s in range(P+1)]
            k0 = d+1
        for k in range(k0, min(d+8, 52)+1):
            if k not in v[dl]: continue
            n = k - k0
            tot = Q(0)
            for s in range(P+1):
                tot = tot + Q(coeffs[s][0], 0, coeffs[s][1], 0).sc(comb(n, s))
            pred = tot*rho0
            cells += 1
            if diag: diag_cells += 1
            got = Q.of(v[dl][k])
            if not (pred - got).iszero():
                bad += 1
                if len(bad_ex) < 5: bad_ex.append((p, q, dl, k))
        if diag:
            # initial level: rho_pp(dl) = rho0 * a_dl(2p, p)
            ar = evE[(dl, 0)](d, p)
            pred = Q(ar[0], 0, ar[1], 0)*rho0
            init_cells += 1
            if not (pred - Q.of(v[dl][d])).iszero(): init_bad += 1
print(f'FULL LAW: cells (all keys d=2..52, delta<=min(6,2d-4), k0..min(d+8,52)) = {cells} (of which diagonal {diag_cells}), bad = {bad}', bad_ex)
print(f'diagonal initial level k=2p vs rho closed form: {init_cells} cells, bad = {init_bad}')
print('time', round(time.time()-t0, 1))
