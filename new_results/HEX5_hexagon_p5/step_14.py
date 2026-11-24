import subprocess
# optimized engine at 50 dps
r1 = subprocess.run(['<project>/.venv/bin/python','step_04.py','50','opt50.txt'],capture_output=True,text=True)
print(r1.stdout[-1400:]); print(r1.stderr[-500:])
