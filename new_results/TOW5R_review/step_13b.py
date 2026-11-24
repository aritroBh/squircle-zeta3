"""step_13b (TOW5R): my own weight-7 kernel check from WTS theory/kern_k.pkl (ground truth, 10-27).
Own Lagrange quintic from k = 8..13; check k = 7 and k = 14 (FRESH level beyond TOW5's range);
own 6th differences on k = 7..13; +-i simple-pole parts empty; k = 6 mismatch count.
Also verify: real-pole residues are real (c=d=0) and +b pole absent, as WTS asserts.
"""
import pickle, os
from fractions import Fraction as Fr
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(6, 15)}
keys6 = [(6,0),(5,1),(4,2),(3,3)]
def lag(ys, xs, x):
    out = Fr(0)
    for a, xa in enumerate(xs):
        La = Fr(1)
        for b, xb in enumerate(xs):
            if b != a: La *= Fr(x-xb, xa-xb)
        out += ys[a]*La
    return out
poles = set()
for key in keys6:
    for k in range(6, 15): poles |= set(K[k][key]['real'].keys())
print('real poles present:', sorted(poles), '(expect 0, a, -a, -b; no +b)', flush=True)
# structurally real?
notreal = sum(1 for key in keys6 for k in range(6,15) for pole in K[k][key]['real']
              for c in K[k][key]['real'][pole][2:] if c != 0)
print('non-real QSI components in real-pole residues:', notreal, flush=True)
n7 = bad7 = n14 = bad14 = 0
for key in keys6:
    pset = set()
    for k in range(7, 14): pset |= set(K[k][key]['real'].keys())
    for pole in pset:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in [8,9,10,11,12,13]]
            n7 += 1
            if lag(ys,[8,9,10,11,12,13],7) != K[7][key]['real'].get(pole,(Fr(0),)*4)[comp]: bad7 += 1
            n14 += 1
            if lag(ys,[8,9,10,11,12,13],14) != K[14][key]['real'].get(pole,(Fr(0),)*4)[comp]: bad14 += 1
print(f'k=7 vs quintic(8..13): {n7} comp, bad {bad7};  k=14 (FRESH) vs quintic: {n14} comp, bad {bad14}', flush=True)
# 6th differences on 7..13
badfd = nfd = 0
for key in keys6:
    pset = set()
    for k in range(7, 14): pset |= set(K[k][key]['real'].keys())
    for pole in pset:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in range(7,14)]
            cur = ys
            for _ in range(6): cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
            nfd += 1
            if any(c != 0 for c in cur): badfd += 1
print(f'6th differences k=7..13: {nfd} comp, bad {badfd}', flush=True)
pmi = sum(1 for key in keys6 for k in range(7,15) if K[k][key]['ipos'] or K[k][key]['ineg'])
print('+-i simple-pole parts nonempty (k=7..14):', pmi, flush=True)
mis6 = tot6 = 0
for key in keys6:
    pset = set()
    for k in list(range(7,14))+[6]: pset |= set(K[k][key]['real'].keys())
    for pole in pset:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr(0),)*4)[comp] for k in [8,9,10,11,12,13]]
            tot6 += 1
            if lag(ys,[8,9,10,11,12,13],6) != K[6][key]['real'].get(pole,(Fr(0),)*4)[comp]: mis6 += 1
print(f'k=6 vs quintic: {tot6} comp, {mis6} mismatched', flush=True)
ok = (bad7 == 0 and bad14 == 0 and badfd == 0 and pmi == 0 and mis6 > 0 and notreal == 0)
print('STEP_13B VERDICT:', 'PASS' if ok else 'FAIL', flush=True)
assert ok
