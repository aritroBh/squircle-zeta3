import subprocess, time, os, sys
cwd = "<project>/research/08_theory/G2_M4alpha"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
PY = sys.executable
t0 = time.time()
r = subprocess.run([PY, G1 + "flat_exact.py", "100", "4"], capture_output=True, text=True, cwd=cwd)
print("time N=100:", time.time() - t0)
print(r.stdout[-400:])
