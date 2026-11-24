import subprocess, time
work = '<project>/research/08_theory/M5a_next_coefficient'
PY = '<project>/.venv/bin/python'
t0 = time.time()
r = subprocess.run(['nice', PY, '-B', 'm3_tails_exact.py', '5'], capture_output=True, text=True, cwd=work)
print("elapsed:", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[:2000])
if 'W3' in r.stdout:
    with open(f'{work}/progress.log', 'a') as f:
        f.write("step: m3_tails_exact.py 5 -> m5_tails_exact.pkl/.out (exact T and W3(5))\n")