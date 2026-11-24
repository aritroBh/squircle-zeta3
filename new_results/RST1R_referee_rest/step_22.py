import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run([PY,'step_02.py'],capture_output=True,text=True)
print("rc:",r.returncode)
print("OUT:",r.stdout)
print("ERR:",r.stderr[-3000:])
