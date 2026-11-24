
import subprocess, time, os, datetime
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
os.chdir(WD)
PY = '<project>/.venv/bin/python'
t0 = time.time()
r = subprocess.run(['nice', '-n', '10', PY, '-B', 'verify_wt6s.py'], capture_output=True, text=True, cwd=WD, timeout=900)
print(f"runtime {time.time()-t0:.0f}s rc={r.returncode}")
print(r.stdout[-1800:])
ok = ('ALL PASS' in r.stdout) and r.returncode == 0
with open('progress.log', 'a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} verify_wt6s.py final run: {'ALL PASS' if ok else 'FAIL'} ({time.time()-t0:.0f}s)\n")
print("VERIFY:", "ALL PASS" if ok else "FAIL")
