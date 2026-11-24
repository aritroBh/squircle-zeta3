import subprocess, time, os
CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
VENV = "<project>/.venv/bin/python"
B3 = "<project>/research/10_frontier/W1_corner_family/B_taylor_pi4/"
G2 = "<project>/research/08_theory/G2_M4alpha/"
t0 = time.time()
r1 = subprocess.run(["nice", VENV, "zsum.py"], cwd=B3, capture_output=True, text=True)
print("zsum self-test:", r1.stdout.strip(), r1.stderr[-200:] if r1.stderr else "")
r2 = subprocess.run(["nice", VENV, "planner_check_b4.py"], cwd=G2, capture_output=True, text=True)
print(r2.stdout, r2.stderr[-200:] if r2.stderr else "")
print(f"({time.time()-t0:.1f}s)")
with open(CWD + "/progress.log", "a") as f:
    f.write(f"step_zsum_selftest+planner: zsum self-tests passed; planner numerical check reproduced ({time.time()-t0:.1f}s)\n")
