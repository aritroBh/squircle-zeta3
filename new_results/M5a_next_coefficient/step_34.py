import subprocess, time
work = '<project>/research/08_theory/M5a_next_coefficient'
PY = '<project>/.venv/bin/python'
t0 = time.time()
r = subprocess.run(['nice', PY, '-B', 'w6_symbol.py'], capture_output=True, text=True, cwd=work)
print("elapsed:", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[:3000])