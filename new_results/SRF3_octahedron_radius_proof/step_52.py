
import subprocess, time
PY = "<project>/.venv/bin/python"
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
t0=time.time()
r = subprocess.run(["nice","-n","10",PY,"verify_srf3.py"],capture_output=True,text=True,cwd=wd,timeout=400)
print(r.stdout); print(r.stderr[-2000:]); print("wall: %.1fs"%(time.time()-t0))
