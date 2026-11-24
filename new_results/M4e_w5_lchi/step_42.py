import subprocess, time
M4E = '<project>/research/08_theory/M4e_w5_lchi'
PY = '<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice', PY, 'nice_lchi2.py', '1200'], cwd=M4E, capture_output=True, text=True, timeout=2700)
print(r.stdout); print("STDERR:", r.stderr[-1500:]); print(f"elapsed {time.monotonic()-t0:.0f}s")
open(M4E+'/progress.log','a').write("step_41: nice_lchi2.py 1200 (saturated kernel)\n")
