import os, subprocess, time
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
with open("progress.log","a") as f:
    f.write("m=3: g1 = ln2 + (1/2)log(1+e^{-4 alpha |s|}) PROVED + verified 1e-12; c3_lay formula validated; [int g2]_zeta = -sqrt3 z3/8 -> c3_lay-zeta = C3/36 (matches paper C3)\n")
PY = "<project>/.venv/bin/python"
t0=time.time()
r = subprocess.run([PY,"av5b_profiles.py","4"],capture_output=True,text=True,timeout=1400)
print(r.stdout, r.stderr[-300:])
print("elapsed",time.time()-t0)
