import subprocess, time, os
CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
VENV = "<project>/.venv/bin/python"
results = {}
for tag in ["mut_a", "mut_b", "mut_c", "mut_d"]:
    t0 = time.time()
    r = subprocess.run(["nice", VENV, "step_11.py"], cwd=os.path.join(CWD, tag),
                       capture_output=True, text=True)
    dt = time.time() - t0
    open(os.path.join(CWD, tag + ".out"), "w").write(r.stdout + "\n--- STDERR ---\n" + r.stderr)
    last = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else "(no output)"
    results[tag] = (last, r.returncode, dt)
    print(f"===== {tag} (rc={r.returncode}, {dt:.1f}s) =====")
    # print only the failing parts + verdict
    for line in r.stdout.splitlines():
        if "False" in line or "residual" in line or "PROVED" in line or "part " in line:
            print(line)
    if r.stderr.strip():
        print("STDERR:", r.stderr[-500:])
with open(CWD + "/progress.log", "a") as fh:
    for tag, (last, rc, dt) in results.items():
        fh.write(f"step_run_{tag}: rc={rc}, {dt:.1f}s, verdict: {last}\n")
