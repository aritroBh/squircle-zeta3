import os, subprocess, time
CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
G2 = "<project>/research/08_theory/G2_M4alpha/"
VENV = "<project>/.venv/bin/python"
RB = os.path.join(CWD, "run_base")
for f in os.listdir(G2):
    if f.startswith("values4_") and f.endswith(".txt"):
        dst = os.path.join(RB, f)
        if not os.path.exists(dst):
            os.symlink(os.path.join(G2, f), dst)
print("symlinked:", [f for f in os.listdir(RB) if f.startswith("values4_")])
t0 = time.time()
r = subprocess.run(["nice", VENV, "step_10.py"], cwd=RB, capture_output=True, text=True)
dt = time.time() - t0
open(CWD + "/run_base_step10.out", "w").write(r.stdout + "\n--- STDERR ---\n" + r.stderr)
print(r.stdout)
print("STDERR:", r.stderr[-300:] if r.stderr else "(empty)")
print(f"rc={r.returncode}, {dt:.1f}s")
with open(CWD + "/progress.log", "a") as fh:
    fh.write("step_run_step10_note: first run rc=1 only because values4_*.txt (relative paths) were absent from run_base; "
             "numerics had already been corroborated by planner_check_b4.py (8.67e-131 at the square)\n")
    fh.write(f"step_run_step10_full: with values4 symlinks, full step_10 (EXACT-1280 + 6-angle numerics) rc={r.returncode}, {dt:.1f}s\n")
