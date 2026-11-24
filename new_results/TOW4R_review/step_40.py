
import os, subprocess, time
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
tc = work + "/tow4_copy"
PY = base + "/.venv/bin/python"
t0 = time.time()
res = {}
for (dl, s) in [(7,0),(7,1),(7,2),(7,3),(8,0),(8,1),(8,2),(8,3),(8,4)]:
    ts = time.time()
    r = subprocess.run(["nice","-n","10",PY,"-B","prove_adj78_fast.py",str(dl),str(s)],
                       capture_output=True, text=True, cwd=tc, timeout=900)
    ok = r.returncode == 0 and "== 0 : True" in r.stdout
    res[(dl,s)] = (ok, round(time.time()-ts,1))
    print(f"ADJ ({dl},{s}): {'PASS' if ok else 'FAIL'} ({round(time.time()-ts,1)}s)", flush=True)
for dl in [7,8]:
    ts = time.time()
    r = subprocess.run(["nice","-n","10",PY,"-B","prove_E3b78.py",str(dl)],
                       capture_output=True, text=True, cwd=tc, timeout=900)
    ok = r.returncode == 0 and "== 0 : True" in r.stdout
    res[('E3',dl)] = (ok, round(time.time()-ts,1))
    print(f"E3 ({dl}): {'PASS' if ok else 'FAIL'} ({round(time.time()-ts,1)}s)", flush=True)
print("ALL:", all(v[0] for v in res.values()), " total", round(time.time()-t0,1))
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} reran 9 ADJ + 2 E3 certificates in tow4_copy: all PASS ({round(time.time()-t0,1)}s)\n")
