
import os, subprocess, time
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
PY='<home>)/.venv/bin/python'
t0 = time.time()
r = subprocess.run([PY, 'planner_sweep_ratio.py', '1.5'], cwd=FIX, capture_output=True, text=True, timeout=7200)
el = time.time()-t0
print(f"({el:.1f}s)")
print(r.stdout)
err = r.stderr[-1500:] if r.stderr else ""
print(err)
lines = r.stdout.strip().split(chr(10))
with open(os.path.join(WORK,'progress.log'),'a') as f:
    f.write("D13s turn8: sweep ratio 1.5 rerun with fixed step_38 (%.0fs): %s | %s | %s\n" % (el, lines[-4] if len(lines)>3 else '', lines[-3] if len(lines)>2 else '', lines[-2] if len(lines)>1 else ''))
