import subprocess, time
work = '<project>/research/08_theory/M5a_next_coefficient'
PY = '<project>/.venv/bin/python'
t0 = time.time()
r = subprocess.run(['nice', PY, '-B', 'step_quad5.py'], capture_output=True, text=True, cwd=work)
print("elapsed:", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[:2000])
if 'M5    =' in r.stdout:
    with open(f'{work}/progress.log', 'a') as f:
        f.write("step: step_quad5.py quadrature of all k=5 pieces + M5 at 100 dps -> quad5_100.json\n")