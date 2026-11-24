code = '''"""step_17 (TOW5R): independent tower-law check from WTS ground truth (tower_k.pkl, k <= 14).
For every key of degree n = 7..12, every depth dl <= min(10, 2n-4), every level k in [k0, 14]
(k0 = n off-diagonal, n+1 diagonal), compare the recorded tower coefficient against
  off-diag: rho(0) * sum_s binom(k-n, s) F(dl,s)(n,p)          (F = tow5lib.forced, PROVED family)
  diagonal: rho_pp(0) * sum_s binom(k-2p-1, s) ehat(dl,s)(p)   (ehat from prove_E3b910)
plus diagonal initial level k = 2p vs a_dl(2p,p)*rho_pp(0).
All arithmetic my own exact Fraction Q(sqrt2,i) (check_data.Q). rho(0) from WTS binteg (towlib.rho).
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE + '/tow5_copy')
_WTS = '<project>/research/08_theory/WTS_weight_k_structure/theory'
sys.path.insert(0, _WTS)
sys.path.insert(0, '<project>/research/08_theory/TOW2_tower_all_depths')
from math import comb
from fractions import Fraction as Fr
import sympy as sp
import tow5lib as T5
import check_data as CD
import towlib
Q = CD.Q
t0 = time.time()
d_, m_ = T5.d, T5.m
# ehat evaluator (import the prove_E3b910 module; it builds ehat on demand)
import importlib.util
spec = importlib.util.spec_from_file_location('PE', _HERE + '/tow5_copy/prove_E3b910.py')
PE = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PE)
p_ = PE.p
evE = {}; evH = {}
for dl in range(0, 11):
    for s in range(0, dl//2+1):
        evE[(dl, s)] = CD.make_eval(T5.forced(dl, s), (d_, m_))
        evH[(dl, s)] = CD.make_eval(sp.cancel(sp.together(PE.ehat(dl, s))), (p_,))
print('evaluators built', round(time.time()-t0,1), flush=True)
def u_wts(i, j, k, delta):
    dd = i+j; r = 2*dd-2-delta
    if r < 2: raise ValueError
    H = towlib.tower(k).get((i,j))
    if H is None: return None
    return H.get(r)
cells = bad = diagc = init = initbad = 0
badex = []
for n in range(7, 13):
    for key in towlib.keys_of_degree(n):
        i0, j0 = key; diag = (i0 == j0)
        rho0 = Q.of((lambda v: (v.p0, v.p1, v.p2, v.p3))(towlib.rho(key, 0)))
        for dl in range(0, min(10, 2*n-4)+1):
            P = dl//2
            if not diag:
                coeffs = [evE[(dl, s)](n, i0) for s in range(P+1)]; k0 = n
            else:
                coeffs = [evH[(dl, s)](i0) for s in range(P+1)]; k0 = n+1
            for k in range(k0, 15):
                v = u_wts(i0, j0, k, dl)
                if v is None: continue
                tot = Q(0)
                for s in range(P+1):
                    tot = tot + Q(coeffs[s][0], 0, coeffs[s][1], 0).sc(comb(k-k0, s))
                pred = tot*rho0
                cells += 1
                if diag: diagc += 1
                got = Q.of((v.p0, v.p1, v.p2, v.p3))
                if not (pred-got).iszero():
                    bad += 1
                    if len(badex) < 5: badex.append((key, dl, k))
            if diag:
                v = u_wts(i0, j0, n, dl)
                if v is not None:
                    ar = evE[(dl, 0)](n, i0)
                    pred = Q(ar[0], 0, ar[1], 0)*rho0
                    init += 1
                    if not (pred - Q.of((v.p0, v.p1, v.p2, v.p3))).iszero(): initbad += 1
    print(f'degree {n} done ({round(time.time()-t0,1)}s)', flush=True)
print(f'STEP_17 VERDICT: WTS-ground-truth law at dl<=10: {cells} cells (diag {diagc}), bad = {bad} {badex}; diag-init {init}, bad {initbad}', flush=True)
assert bad == 0 and initbad == 0
'''
wd = '<project>/research/08_theory/TOW5R_review'
open(wd+'/step_17.py','w').write(code)
import subprocess, datetime
PY = '<project>/.venv/bin/python'
t0 = datetime.datetime.now()
r = subprocess.run(['nice','-n','10',PY,'step_17.py'], cwd=wd, capture_output=True, text=True, timeout=1700)
print(r.stdout[-1200:])
print('STDERR tail:', r.stderr[-400:], 'rc:', r.returncode, 'elapsed:', datetime.datetime.now()-t0)
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_17: independent WTS-ground-truth law dl<=10, deg 7..12, k<=14 rc={r.returncode} ({(datetime.datetime.now()-t0).seconds}s)\n")
