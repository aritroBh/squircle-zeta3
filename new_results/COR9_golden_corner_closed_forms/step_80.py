import subprocess, time
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
t0 = time.time()
proc = subprocess.run(['<project>/.venv/bin/python', WD+'verify_cor9.py'],
                      capture_output=True, text=True, timeout=280)
print(proc.stdout[-4000:])
print("STDERR:", proc.stderr[-1500:] if proc.stderr else "(none)")
print(f"wall time: {time.time()-t0:.0f}s, returncode {proc.returncode}")
