
import os, subprocess, time
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
tc = work + "/tow4_copy"
PY = base + "/.venv/bin/python"
t0=time.time()
for (dl,s) in [(7,3),(8,4),(8,3),(7,2)]:
    ts=time.time()
    try:
        r = subprocess.run(["nice","-n","10",PY,"-B","prove_adj78.py",str(dl),str(s)],
                           capture_output=True, text=True, cwd=tc, timeout=1500)
        out=(r.stdout.strip().splitlines() or [""])[-1]
        print(f"SLOW ({dl},{s}): rc={r.returncode} '{out[:95]}' ({round(time.time()-ts,1)}s)", flush=True)
    except subprocess.TimeoutExpired:
        print(f"SLOW ({dl},{s}): TIMEOUT >1500s", flush=True)
        break
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} fast-vs-slow cross-validation prove_adj78.py (7,3),(8,4),(8,3),(7,2)\n")
