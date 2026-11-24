import os, subprocess, time
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
PY = "<project>/.venv/bin/python"
t0=time.time()
r = subprocess.run([PY,"av5b_profiles.py","3"],capture_output=True,text=True,timeout=1400)
print(r.stdout, r.stderr[-400:])
print("elapsed",time.time()-t0)
