import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_02.py','150','p5_run150.txt'],capture_output=True,text=True)
print(r.stdout); print(r.stderr)
