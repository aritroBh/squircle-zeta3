import subprocess, time, os
CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
VENV = "<project>/.venv/bin/python"
t0 = time.time()
r = subprocess.run(["nice", VENV, "step_11.py"], cwd=os.path.join(CWD, "run_base"),
                   capture_output=True, text=True)
dt = time.time() - t0
open(os.path.join(CWD, "run_base.out"), "w").write(r.stdout + "\n--- STDERR ---\n" + r.stderr)
print(r.stdout)
print("stderr tail:", r.stderr[-500:] if r.stderr else "(empty)")
print(f"returncode={r.returncode}, wall time {dt:.1f}s")
with open(os.path.join(CWD, "progress.log"), "a") as f:
    f.write(f"step_run_base: pristine step_11.py run in run_base/, rc={r.returncode}, {dt:.1f}s, verdict line: {r.stdout.strip().splitlines()[-1]}\n")
