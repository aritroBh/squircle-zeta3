import subprocess, time
M4E = '<project>/research/08_theory/M4e_w5_lchi'
PY = '<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice', PY, 'verify_nice.py', '1500'], cwd=M4E, capture_output=True, text=True, timeout=2700)
# print only the PASS/value lines
for line in r.stdout.splitlines():
    if line.startswith(('1.','2a','2b','3.','4.','ALL')) or 'PASS' in line:
        print(line)
print("STDERR:", r.stderr[-800:]); print(f"elapsed {time.monotonic()-t0:.0f}s")
open(M4E+'/progress.log','a').write("step_49: verify_nice.py 1500\n")
