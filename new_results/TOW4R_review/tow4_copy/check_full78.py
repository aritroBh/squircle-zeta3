"""check_full78.py -- full tower-law exact confirmation at depths <= 8 vs TOWN1 towers.
Every data cell (key, delta <= min(8,2d-4), k in [k0, min(d+8,52)]), k0 = d (off-diagonal) /
d+1 (diagonal stable), compared with the Newton-form prediction from the PROVED e-family:
off-diagonal tow4lib.forced; diagonal ehat (prove_E3b78). Plus diagonal initial level k=2p
vs rho closed form a_delta(2p,p). Mirrors TOW3 referee check_full.py (38939 cells at delta<=6).
"""
import sys, pickle, time
from math import comb
sys.path.insert(0, '<project>/research/08_theory/TOW4R_review/tow4_copy')
import sympy as sp
import check_data as CD
import tow4lib as T4
import prove_E3b78 as PE
Q = CD.Q
t0 = time.time()
H = pickle.load(open(CD.PKL, 'rb')); H.pop('_meta')
d_, m_ = T4.d, T4.m
p_ = PE.p
evE = {}; evH = {}
for dl in range(0, 9):
    for s in range(0, dl//2+1):
        evE[(dl, s)] = CD.make_eval(T4.forced(dl, s), (d_, m_))
        evH[(dl, s)] = CD.make_eval(sp.cancel(sp.together(PE.ehat(dl, s))), (p_,))
print('evaluators built', round(time.time()-t0, 1), flush=True)
cells = 0; bad = 0; diag_cells = 0; init_cells = 0; init_bad = 0
bad_ex = []
for (p, q), v in sorted(H.items()):
    d = p+q
    diag = (p == q)
    for dl in range(0, min(8, 2*d-4)+1):
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
            ar = evE[(dl, 0)](d, p)
            pred = Q(ar[0], 0, ar[1], 0)*rho0
            init_cells += 1
            if not (pred - Q.of(v[dl][d])).iszero(): init_bad += 1
print(f'FULL LAW delta<=8: cells = {cells} (diagonal {diag_cells}), bad = {bad}', bad_ex, flush=True)
print(f'diagonal initial k=2p vs a_delta(2p,p): {init_cells} cells, bad = {init_bad}', flush=True)
print('time', round(time.time()-t0, 1), flush=True)
