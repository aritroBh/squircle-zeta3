import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice',PY,'verify_rul1.py'], capture_output=True, text=True, timeout=1800)
print(r.stdout[-5000:]); print("STDERR:", r.stderr[-1500:])
