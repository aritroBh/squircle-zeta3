import os, subprocess, time
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
PY = "<project>/.venv/bin/python"
t0 = time.time()
r = subprocess.run([PY, "verify_av5b.py"], capture_output=True, text=True, timeout=890)
print(r.stdout)
print("ERR:", r.stderr[-800:])
print("elapsed", time.time()-t0)
