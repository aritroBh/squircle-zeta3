
import subprocess, os, time
os.chdir("research/08_theory/D13e_certify")
PY = "<home>)/.venv/bin/python"
r = subprocess.run([PY, "dbg1.py"], capture_output=True, text=True, timeout=300)
print("STDOUT:", r.stdout)
print("STDERR:", r.stderr[:2000])
print("="*30, "step_35.py", "="*30)
print(open("step_35.py").read())
