import subprocess, sys, os, time
cwd = "<project>/research/08_theory/G2_M4alpha"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
PY = sys.executable
# background run of the reference implementation
log = open(os.path.join(cwd, "flat_exact_N400_K4.out"), "w")
p = subprocess.Popen([PY, G1 + "flat_exact.py", "400", "4"], stdout=log, stderr=subprocess.STDOUT, cwd=cwd)
print("launched pid", p.pid)
# save pid
open(os.path.join(cwd, ".flat400.pid"), "w").write(str(p.pid))
