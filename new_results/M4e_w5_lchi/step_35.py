import subprocess, time
M4E = '<project>/research/08_theory/M4e_w5_lchi'
PY = '<project>/.venv/bin/python'
t0 = time.monotonic()
r = subprocess.run(['nice', PY, 'nice_lchi.py', '1200'], cwd=M4E,
                   capture_output=True, text=True, timeout=2400)
print(r.stdout)
print("STDERR:", r.stderr[-3000:])
print(f"elapsed {time.monotonic()-t0:.0f}s")
open(M4E+'/progress.log','a').write("step_33: nice_lchi.py 1200 run\n")
