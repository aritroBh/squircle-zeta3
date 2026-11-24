import os, subprocess, time
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
PY = "<project>/.venv/bin/python"
t0=time.time()
r = subprocess.run([PY,"av5b_newfam.py","run","0.048,0.024,0.012,0.006,0.003,0.0015,0.00075,0.000375,0.0001875"],
                   capture_output=True,text=True,timeout=1400)
print(r.stdout, r.stderr[-300:])
print("elapsed",time.time()-t0)
