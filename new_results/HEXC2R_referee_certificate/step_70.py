
import subprocess
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "steps/step_09.py"], capture_output=True, text=True, timeout=3000)
print(r.stdout)
print(r.stderr[-1000:] if r.returncode else "OK")
