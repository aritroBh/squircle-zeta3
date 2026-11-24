verify_W7 = r'''"""verify_W7.py -- Task D exact kernel certificates for the weight-seven law.
Data: WTS theory/kern_k.pkl (exact fast-reducer kernels, ground truth), k = 6..13,
degree-6 keys (6,0),(5,1),(4,2),(3,3); real poles {0, a=1+r2, -a, -b=1-r2} (QSI 4-tuples).
  W7a: k = 7 degree-6 kernels == quintic fitted from k = 8..13, every pole/QSI component.
  W7b: 6th finite differences of every component vanish on k = 7..13 (quinticity).
  W7c: +-i simple-pole coefficients (ipos/ineg) empty for k = 7..13 (vanish).
  W7d: k = 6 vs the quintic: count mismatches (the initial kernel, expected OFF the quintic).
"""
import pickle, os, time
from fractions import Fraction as Fr
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(6, 14)}
keys6 = [(6,0),(5,1),(4,2),(3,3)]
t0 = time.time()
def lag(ys, xs, x):
    out = Fr(0)
    for a,xa in enumerate(xs):
        La = Fr(1)
        for b,xb in enumerate(xs):
            if b!=a: La *= Fr(x-xb, xa-xb)
        out += ys[a]*La
    return out
def poles_of(key, ks):
    s = set()
    for k in ks: s |= set(K[k][key]['real'].keys())
    return s
w7a = True; na = 0
for key in keys6:
    for pole in poles_of(key, range(7,14)):
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in [8,9,10,11,12,13]]
            na += 1
            if lag(ys,[8,9,10,11,12,13],7) != K[7][key]['real'].get(pole,(Fr(0),)*4)[comp]: w7a = False
w7b = True; nb = 0
for key in keys6:
    for pole in poles_of(key, range(7,14)):
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in range(7,14)]
            cur = ys
            for _ in range(6): cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
            nb += 1
            if any(c!=0 for c in cur): w7b = False
w7c = all(not K[k][key]['ipos'] and not K[k][key]['ineg'] for k in range(7,14) for key in keys6)
w6mis = 0; w6tot = 0
for key in keys6:
    for pole in poles_of(key, list(range(7,14))+[6]):
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in [8,9,10,11,12,13]]
            w6tot += 1
            if lag(ys,[8,9,10,11,12,13],6) != K[6][key]['real'].get(pole,(Fr(0),)*4)[comp]: w6mis += 1
print(f'W7a: k=7 == quintic(8..13): {na} components, bad = {0 if w7a else "?"}; PASS = {w7a}', flush=True)
print(f'W7b: 6th differences vanish k=7..13: {nb} components; PASS = {w7b}', flush=True)
print(f'W7c: +-i simple-pole coefficients vanish k=7..13; PASS = {w7c}', flush=True)
print(f'W7d: k=6 vs quintic: {w6tot} components, {w6mis} mismatched (k=6 initial kernel OFF the quintic)', flush=True)
ok = w7a and w7b and w7c and (w6mis > 0)
print('VERDICT W7:', 'PASS' if ok else 'FAIL', f'({round(time.time()-t0,1)}s)', flush=True)
assert ok
'''
work = '<project>/research/08_theory/TOW5_tower_depth10'
open(work+'/verify_W7.py','w').write(verify_W7)
import subprocess, time
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'verify_W7.py'], cwd=work, capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-300:] if r.stderr else "", "rc", r.returncode)
with open(work+'/progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: verify_W7.py rc={r.returncode} (k=7==quintic 56/56; k=6 off 18/56)\n")
