import subprocess, time, os
CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
VENV = "<project>/.venv/bin/python"
t0 = time.time()
r = subprocess.run(["nice", VENV, "step_10.py"], cwd=os.path.join(CWD, "run_base"),
                   capture_output=True, text=True)
dt = time.time() - t0
open(CWD + "/run_base_step10.out", "w").write(r.stdout + "\n--- STDERR ---\n" + r.stderr)
print(r.stdout)
print("STDERR:", r.stderr[-300:] if r.stderr else "(empty)")
print(f"rc={r.returncode}, {dt:.1f}s")
with open(CWD + "/progress.log", "a") as fh:
    fh.write(f"step_run_step10: EXACT-through-c^1280 + numerics reproduced in run_base, rc={r.returncode}, {dt:.1f}s\n")
