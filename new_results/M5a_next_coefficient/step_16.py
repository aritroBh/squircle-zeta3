import subprocess, time, os
base = '<project>'
work = f'{base}/research/08_theory/M5a_next_coefficient'
PY = f'{base}/.venv/bin/python'
t0 = time.time()
r = subprocess.run(['nice', PY, '-B', 'ibp_reduce.py', '5'], capture_output=True, text=True, cwd=work)
print("elapsed:", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[:2000])
open(f'{work}/ibp_reduce_k5.out', 'w').write(r.stdout)
with open(f'{work}/progress.log', 'a') as f:
    f.write(f"step: ran ibp_reduce.py 5 -> ibp_reduce_k5.out ({time.time()-t0:.0f}s)\n")