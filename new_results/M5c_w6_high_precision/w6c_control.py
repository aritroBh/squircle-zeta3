"""M5c full-lattice positive control: plant a known target relation in the 179-column W6 lattice
and check LLL recovers it. Target := sum of known small-coefficient combos of kv/sqrt2-kv columns.
Usage: -B w6c_control.py <digits>
"""
import json, math, pickle, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import mpmath as mp
from w6c_lll import build, col_values, set_prec
from m5b_common import log

def main():
    digits = int(sys.argv[1]) if len(sys.argv) > 1 else 1300
    d = build('')
    nk = len(d['kv'])
    kept = list(range(2*nk))
    set_prec(digits+120)
    vals = col_values(d, kept)   # [target, kv0, r2 kv0, kv1, r2 kv1, ...]
    # planted relation: target = 3*kv_2 - 5*sqrt2 kv_2 + 7*kv_5 + 11*sqrt2 kv_5 - 13*kv_40 + 2*sqrt2 kv_40
    plant = {2: 3, 3: -5, 10: 7, 11: 11, 80: -13, 81: 2}
    planted_target = mp.fsum(mp.mpf(c)*vals[1+j] for j, c in plant.items())
    vals = [planted_target] + vals[1:]
    n = len(vals)
    C = mp.mpf(10)**digits
    rows = [[int(i == j) for j in range(n)]+[int(mp.nint(C*v))] for i, v in enumerate(vals)]
    log(f'control LLL start: rows={n}, digits={digits}, planted={plant}')
    t0 = time.monotonic()
    red = flint.fmpz_mat(rows).lll(gram='exact')
    dt = time.monotonic()-t0
    expected = [1]+[-plant.get(j, 0) for j in range(2*nk)]
    found = None
    set_prec(digits+200)
    for i in range(n):
        vec = [int(red[i, j]) for j in range(n)]
        if vec == expected or vec == [-x for x in expected]:
            found = vec
            break
    res = mp.fsum(mp.mpf(c)*v for c, v in zip(expected, vals))
    log(f'control LLL done: seconds={dt:.1f}; planted relation recovered={found is not None}; residual={mp.nstr(res, 6)}')
    if found is None:
        # show top rows for debugging
        diag = []
        for i in range(n):
            vec = [int(red[i, j]) for j in range(n)]
            r = mp.fsum(mp.mpf(c)*v for c, v in zip(vec, vals))
            diag.append((abs(r), i, max(map(abs, vec))))
        diag.sort()
        for r, i, mx in diag[:5]:
            log(f'  control row {i}: maxcoef={mx} residual={mp.nstr(r, 6)}')
    assert found is not None, 'FULL-LATTICE CONTROL FAILED'
    log('full-lattice planted control: PASS')

if __name__ == '__main__':
    main()
