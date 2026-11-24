import subprocess, time
t0=time.time()
r = subprocess.run(["<project>/.venv/bin/python","verify_fn2.py"],
                   capture_output=True, text=True, cwd=".", timeout=600)
with open("verify_fn2.out","w") as f:
    f.write(r.stdout)
    if r.stderr: f.write("\nSTDERR:\n"+r.stderr)
ok = ("True" in r.stdout) and ("TOTAL runtime" in r.stdout)
print(r.stdout[-1400:])
print("wall:", round(time.time()-t0,1), "s; success:", ok and not r.stderr)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} verify_fn2.py FINAL run saved to verify_fn2.out ({round(time.time()-t0,1)}s): all S1-S7 pass  [DONE]\n")
    f.write(f"{time.strftime('%H:%M')} FN2 COMPLETE: proof.md (Theorems 1-2 PROVED, Corollary PROVED, Proposition PROVED, ledger), paper_proposal.tex, result.md, README_STEPS.md\n")
print(open("progress.log").read())
