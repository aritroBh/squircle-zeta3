
import subprocess
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "steps/step_10.py"], capture_output=True, text=True, timeout=1800)
print(r.stdout)
print(r.stderr[-3000:] if r.returncode else "OK")
