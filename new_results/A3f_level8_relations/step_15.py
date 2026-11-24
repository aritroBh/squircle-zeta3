code = r'''"""step_03_numeric_w123.py: evaluate every audited relation row (weights 1-3) at ~60 digits.
acb rigorous balls, ctx.prec=200 (~60 digits); assert each |value| < 1e-40 (midpoint)."""
import sys, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
import check_conversion as CC
from check_conversion import evaluate, approx, mp

def main():
    for n in (1, 2, 3):
        with open(f'audit_weight{n}.pickle', 'rb') as f:
            data = pickle.load(f)
        rows = data['rows']
        t0 = time.monotonic()
        worst = mp.mpf(0); bad = 0
        for row in rows:
            v = abs(approx(evaluate(dict(row))))
            worst = max(worst, v)
            if v > mp.mpf('1e-40'):
                bad += 1
                if bad <= 3: print('BAD', n, v, row[:5])
        print(f'weight {n}: rows={len(rows)} bad={bad} worst={mp.nstr(worst,8)} ({time.monotonic()-t0:.1f}s)', flush=True)
        assert bad == 0
    print('NUMERIC AUDIT w1-3 PASS at ~60 digits')

main()
'''
open('step_03_numeric_w123.py','w').write(code)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-B','step_03_numeric_w123.py'],
                   capture_output=True, text=True, timeout=2400)
print(r.stdout); print('STDERR:', r.stderr[-2000:])
