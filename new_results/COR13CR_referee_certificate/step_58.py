import subprocess, time, os
T0 = time.monotonic()
py = '<project>/.venv/bin/python'
r = subprocess.run(['nice', py, 'verify_cor13cr.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout)
print("STDERR:", r.stderr[-800:])
print("rc:", r.returncode, f"({time.monotonic()-T0:.0f}s)")
