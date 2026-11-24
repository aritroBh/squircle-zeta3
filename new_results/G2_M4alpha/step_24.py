import subprocess, time, os
cwd = "<project>/research/08_theory/G2_M4alpha"
os.chdir(cwd)
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
t0 = time.time()
r = subprocess.run(["python3", G1 + "flat_exact.py", "30", "4"], capture_output=True, text=True, cwd=cwd)
print("time", time.time() - t0)
print(r.stdout[-3000:])
print(r.stderr[-1000:])
