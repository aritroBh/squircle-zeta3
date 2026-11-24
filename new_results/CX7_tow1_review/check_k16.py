"""Independent fresh integrand/reducer check; run with the project Python -B."""
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'WT2_weight2_all_k'))
import fastred as f
from buildintegrand import build_integrand_pf
from fractions import Fraction as F

start = time.monotonic()
k = 16
initial = build_integrand_pf(k)
print(f'k={k}: fresh integrand built in {time.monotonic()-start:.2f}s', flush=True)
reduced = f.reduce_k_fast(k, integrand_pf=initial)
expected20, expected11 = f.pf_zero(), f.pf_zero()
expected20[0] = {f.ZERO: f.QSI(F(1,4), F(-1,4)), f.AQ: f.QSI(F(1,4)),
                 -f.AQ: f.QSI(F(2*k-7,4)), -f.BQ: f.QSI(F(-(k-2),2))}
expected11[0] = {f.ZERO: f.QSI(F(2*k-5,4), F(1,4)), -f.AQ: f.QSI(F(-(2*k-5),2))}
for key, expected in (((2,0), expected20), ((1,1), expected11)):
    actual = reduced.get(key, f.pf_zero())
    diff = f.pf_sub(actual, expected)
    assert not diff[0] and all(x.is_zero() for part in diff[1:3] for x in part) and diff[3].is_zero(), key
    assert not actual[1] and not actual[2] and actual[3].is_zero(), key
    print(f'PASS k={k} key={key}: exact zero residual; no imaginary poles or polynomial part', flush=True)
    print(f'  kernel = {f.pf_to_sympy(actual)}', flush=True)
print(f'PASS: both fresh k={k} kernels; total {time.monotonic()-start:.2f}s', flush=True)
