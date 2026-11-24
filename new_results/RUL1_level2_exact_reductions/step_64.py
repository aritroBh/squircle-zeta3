import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice',PY,'rul1_step04.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout); print(r.stderr[-800:])
