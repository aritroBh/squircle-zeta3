
import subprocess, time, os, json
PY = "<project>/.venv/bin/python"
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
t0=time.time()
r = subprocess.run(["nice","-n","10",PY,"verify_srf3.py"],capture_output=True,text=True,cwd=wd,timeout=400)
print(r.stdout)
print("wall: %.1fs"%(time.time()-t0))
# confirm K=120 artifact intact
h = json.load(open(wd+"/SRF3_hi_coeffs_K120.json"))
print("K120 artifact K =", h['K'], "len eta:", len(h['eta']))
