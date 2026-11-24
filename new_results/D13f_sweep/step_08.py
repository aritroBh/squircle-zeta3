
import subprocess, os, time
os.chdir("research/08_theory/D13f_sweep")
PY = "<home>)/.venv/bin/python"
t0 = time.time()
r = subprocess.run([PY, "step_34.py", "100", "400"], capture_output=True, text=True)
print("elapsed", time.time()-t0)
print("STDOUT:", r.stdout)
print("STDERR:", r.stderr[:3000])
