import re
src = open('step_05_weight2_values.py').read()
src = src.replace('ctx.cap = 1150', 'ctx.cap = 2048').replace('CC.N = 1150', 'CC.N = 2048')
src = src.replace('series truncation ~ 2^-1150 ~ 1e-346', 'series truncation ~ 0.654^2048 ~ 1e-375 (slowest pole 1-zeta8)')
open('step_05_weight2_values.py','w').write(src)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-B','step_05_weight2_values.py'],
                   capture_output=True, text=True, timeout=2400)
print(r.stdout[-1500:]); print('STDERR:', r.stderr[-1000:])
r2 = subprocess.run(['<home>)/.venv/bin/python','-B','step_06_weight2_lll.py'],
                   capture_output=True, text=True, timeout=2400)
print(r2.stdout); print('STDERR2:', r2.stderr[-1000:])
