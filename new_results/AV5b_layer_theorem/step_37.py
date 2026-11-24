import os, subprocess, time
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
PY = "<project>/.venv/bin/python"
ds = ",".join([mp.mpf('0.048')/2**k for k in range(10)]) if False else "0.048,0.024,0.012,0.006,0.003,0.0015,0.00075,0.000375,0.0001875,0.00009375"
t0=time.time()
for mm in ["2","3"]:
    r = subprocess.run([PY,"step_24.py",mm,"0",ds],capture_output=True,text=True,timeout=1400)
    print(r.stdout)
print("elapsed",time.time()-t0)
