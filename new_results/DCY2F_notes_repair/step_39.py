import subprocess, time, os
CWD = "<project>/research/08_theory/DCY2F_notes_repair"
PY = "<project>/.venv/bin/python"
t0 = time.monotonic()
r = subprocess.run(["nice","-n","10",PY,"verify_dcy2_repaired.py"],
                   cwd=CWD, capture_output=True, text=True, timeout=900)
dt = time.monotonic()-t0
print(r.stdout)
print("STDERR:", r.stderr[-3000:] if r.stderr else "(none)")
print("RETURN CODE:", r.returncode, f" wall={dt:.1f}s")
