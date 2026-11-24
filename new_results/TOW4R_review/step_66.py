
import os, subprocess, time
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
tc = work + "/tow4_copy"
PY = base + "/.venv/bin/python"
t0=time.time()
r = subprocess.run(["nice","-n","10",PY,"-B","verify_tow4.py"], capture_output=True, text=True, cwd=tc, timeout=2400)
print(r.stdout)
print("STDERR:", r.stderr[-300:] if r.stderr else "")
print("elapsed", round(time.time()-t0,1))
open(work+"/verify_tow4_rerun.log","w").write(r.stdout + ("\nSTDERR:\n"+r.stderr if r.stderr else ""))
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} full verify_tow4.py rerun in tow4_copy: {r.stdout.strip().splitlines()[-4].strip()} ({round(time.time()-t0,1)}s)\n")
