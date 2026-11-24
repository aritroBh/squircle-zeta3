import subprocess, sys, os
cwd = "<project>/research/08_theory/G2_M4alpha"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
PY = sys.executable
log = open(os.path.join(cwd, "flat_exact_N520_K4.out"), "w")
p = subprocess.Popen([PY, G1 + "flat_exact.py", "520", "4"], stdout=log, stderr=subprocess.STDOUT, cwd=cwd)
print("launched flat_exact 520 pid", p.pid)
