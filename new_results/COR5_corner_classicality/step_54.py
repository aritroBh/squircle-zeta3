import subprocess, sys
g1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
r = subprocess.run([sys.executable, "m3_hpl_check.py", "60", "2pi_3"], cwd=g1, capture_output=True, text=True, timeout=1200)
print(r.stdout)
print(r.stderr[-1500:] if r.stderr else "")
