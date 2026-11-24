step28 = r'''import sys, time, json
sys.dont_write_bytecode = True
sys.path.insert(0, '../HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from cor13c_common import target2, target3
import exact_linalg as XL
from fractions import Fraction as Fr
t0 = time.monotonic()
lower = {}
cert = {}
for n in (1,2,3):
    rows, _ = cmzv6.make_relations(n, lower)
    lower[n] = rows
    cols = cmzv6.words(n)
    if n < 2: continue
    T = {2: target2, 3: target3}[n]()
    prow, r = XL.select_independent_rows(rows, cols)
    sel = [rows[i] for i in prow]
    pcols = XL.select_pivot_cols(sel, cols)
    y = XL.exact_solve(sel, pcols, cols, T)
    ok, bad = XL.verify_combination(sel, y, cols, T)
    print(f"weight {n}: rank={r} solved; exact verify ok={ok} ({time.monotonic()-t0:.1f}s)")
    assert ok, bad
    hn = [max(len(str(abs(c.numerator))), len(str(c.denominator))) for c in y]
    print(f"   multipliers {len(y)}, max height {max(hn)} digits")
    cert[str(n)] = dict(cols=[list(w) for w in cols], prow=prow,
                        y=[[str(c.numerator), str(c.denominator)] for c in y])
json.dump(cert, open('cert13_w23.json', 'w'))
import os
print('cert13_w23.json', os.path.getsize('cert13_w23.json'), 'bytes')
with open('progress.log','a') as f:
    f.write('step_28: exact certs w2/w3 solved+verified, saved cert13_w23.json\n')
'''
open('step_28.py','w').write(step28)
import subprocess, time
t0 = time.time()
r = subprocess.run(['<project>/.venv/bin/python', 'step_28.py'], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-3000:] if r.stderr else "")
print(f"elapsed {time.time()-t0:.1f}s")
