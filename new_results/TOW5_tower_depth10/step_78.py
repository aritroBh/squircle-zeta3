check_ext = r'''"""check_full910_ext.py -- EXTENDED full tower-law confirmation at depths <= 10.
Data: data_towers_D10_K30_deg18.pkl (levels k=7..30, degrees 7..18, 16291 cells, exact Q(sqrt2,i)).
Same Newton-polynomial check as check_full910.py (forced family off-diagonal, ehat diagonal),
plus diagonal initial level k = 2p vs a_delta(2p,p).
"""
import sys, os, pickle, time
sys.dont_write_bytecode = True
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work)
from math import comb
from fractions import Fraction as Fr
import sympy as sp
import tow5lib as T5
import prove_E3b910 as PE
import check_data as CD
Q = CD.Q
t0 = time.time()
TOW = pickle.load(open(os.path.join(work, 'data_towers_D10_K30_deg18.pkl'), 'rb'))
d_, m_ = T5.d, T5.m
p_ = PE.p
evE = {}; evH = {}
for dl in range(0, 11):
    for s in range(0, dl//2+1):
        evE[(dl, s)] = CD.make_eval(T5.forced(dl, s), (d_, m_))
        evH[(dl, s)] = CD.make_eval(sp.cancel(sp.together(PE.ehat(dl, s))), (p_,))
print('evaluators built', round(time.time()-t0, 1), flush=True)
cells = 0; bad = 0; diag_cells = 0; init_cells = 0; init_bad = 0; bad_ex = []
bykey = {}
for k in TOW:
    for key in TOW[k]:
        for delta, val in TOW[k][key].items():
            bykey.setdefault(key, {}).setdefault(delta, {})[k] = val
for (p, q), vd in sorted(bykey.items()):
    n = p+q; diag = (p == q)
    for dl in range(0, min(10, 2*n-4)+1):
        kd = vd.get(dl, {})
        rho0 = Q.of(vd[0][n])
        P = dl//2
        if not diag:
            coeffs = [evE[(dl, s)](n, p) for s in range(P+1)]
            k0 = n
        else:
            coeffs = [evH[(dl, s)](p) for s in range(P+1)]
            k0 = n+1
        for k in range(k0, 31):
            if k not in kd: continue
            nn = k - k0
            tot = Q(0)
            for s in range(P+1):
                tot = tot + Q(coeffs[s][0], 0, coeffs[s][1], 0).sc(comb(nn, s))
            pred = tot*rho0
            cells += 1
            if diag: diag_cells += 1
            got = Q.of(kd[k])
            if not (pred - got).iszero():
                bad += 1
                if len(bad_ex) < 5: bad_ex.append((p, q, dl, k))
        if diag and n in kd:
            ar = evE[(dl, 0)](n, p)
            pred = Q(ar[0], 0, ar[1], 0)*rho0
            init_cells += 1
            if not (pred - Q.of(kd[n])).iszero(): init_bad += 1
print(f'FULL LAW (extended) delta<=10: cells = {cells} (diagonal {diag_cells}), bad = {bad}', bad_ex, flush=True)
print(f'diagonal initial k=2p vs a_delta(2p,p): {init_cells} cells, bad = {init_bad}', flush=True)
print('time', round(time.time()-t0, 1), flush=True)
assert bad == 0 and init_bad == 0
'''
work = '<project>/research/08_theory/TOW5_tower_depth10'
open(work+'/check_full910_ext.py','w').write(check_ext)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'check_full910_ext.py'], cwd=work, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-500:] if r.stderr else "", "elapsed", round(time.time()-t0,1))
with open(work+'/progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: check_full910_ext.py rc={r.returncode} ({round(time.time()-t0,1)}s)\n")
