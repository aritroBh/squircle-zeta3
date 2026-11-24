import os, subprocess, time
CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
VENV = "<project>/.venv/bin/python"
os.makedirs(CWD + "/flat_run", exist_ok=True)
t0 = time.time()
r = subprocess.run(["nice", VENV, G1 + "flat_exact.py", "40", "4"], cwd=CWD + "/flat_run",
                   capture_output=True, text=True)
print(r.stdout[-1200:])
print("STDERR:", r.stderr[-300:] if r.stderr else "(empty)")
print(f"rc={r.returncode}, {time.time()-t0:.1f}s")
with open(CWD + "/progress.log", "a") as fh:
    fh.write(f"step_flat_exact40: G1 flat_exact.py 40 4 (definition of mu_4,n) run, rc={r.returncode}, {time.time()-t0:.1f}s\n")
