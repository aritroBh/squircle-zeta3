"""Exact k=4 extrapolation, simple poles, and optional fresh full reduction."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import pickle
import json
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'WT2_weight2_all_k'))
import fastred as f

def load(k, key):
    with (HERE.parent / 'WTS_weight_k_structure' / 'theory' / f'kern_{k}.pkl').open('rb') as src:
        p = pickle.load(src)[key]
    return [{f.QSI(*r): f.QSI(*c) for r,c in p['real'].items()},
            [f.QSI(*c) for c in p['ipos']], [f.QSI(*c) for c in p['ineg']], f.QSI(*p['const'])]

def iszero(p):
    return not p[0] and all(c.is_zero() for part in p[1:3] for c in part) and p[3].is_zero()

def polynomial(k, key):
    p5, p6, p7 = [load(j,key) for j in (5,6,7)]
    first = f.pf_sub(p6,p5)
    second = f.pf_add(f.pf_sub(p7,f.pf_scale(p6,2)),p5)
    return f.pf_add(f.pf_add(p5,f.pf_scale(first,k-5)),f.pf_scale(second,F((k-5)*(k-6),2)))

for key in ((3,0),(2,1)):
    assert iszero(f.pf_sub(load(4,key),polynomial(4,key))), key
    for k in (4,5,6,7):
        p = load(k,key)
        assert all(c.is_zero() for part in p[1:3] for c in part) and p[3].is_zero(), (key,k)
    print(f'PASS {key}: k=4 equals extrapolation 3*K5 - 3*K6 + K7; all +/-i coefficients and constants zero at k=4,5,6,7',flush=True)
    p4,p5,p6 = [load(j,key) for j in (4,5,6)]
    first = f.pf_sub(p5,p4)
    second = f.pf_add(f.pf_sub(p6,f.pf_scale(p5,2)),p4)
    print('K4 =', f.pf_to_sympy(p4),flush=True)
    print('Delta K4 =', f.pf_to_sympy(first),flush=True)
    print('Delta^2 K4 =', f.pf_to_sympy(second),flush=True)
    deltas=json.loads((HERE.parent/'WT4_weight4_all_k'/'deltas_deg3.json').read_text())
    poles={'0':f.ZERO,'a':f.AQ,'-a':-f.AQ,'b':f.BQ,'-b':-f.BQ}
    for label,actual in [('D1',first),('D2',second)]:
        expected=f.pf_zero()
        expected[0]={poles[nm]:f.QSI(F(c[0]),F(c[1])) for nm,c in deltas[f'{key[0]},{key[1]}'][label]}
        assert iszero(f.pf_sub(actual,expected)),(key,label)
    print(f'PASS {key}: first/second difference kernels equal the WT4 D1/D2 certificate targets',flush=True)

if len(sys.argv)>1:
    from buildintegrand import build_integrand_pf
    k = int(sys.argv[1])
    start = time.monotonic()
    initial = build_integrand_pf(k)
    print(f'Fresh k={k} integrand built in {time.monotonic()-start:.2f}s',flush=True)
    reduced = f.reduce_k_fast(k,integrand_pf=initial)
    for key in ((3,0),(2,1)):
        actual = reduced.get(key,f.pf_zero())
        assert iszero(f.pf_sub(actual,polynomial(k,key))), (k,key)
        assert all(c.is_zero() for part in actual[1:3] for c in part) and actual[3].is_zero()
        print(f'PASS fresh k={k} {key}: exact kernel residual zero; no +/-i poles or polynomial part',flush=True)
    print(f'PASS fresh reduction, total {time.monotonic()-start:.2f}s',flush=True)
