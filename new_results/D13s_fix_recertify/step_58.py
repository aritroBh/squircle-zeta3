
import subprocess, os, json
ROOT = 'research/08_theory/D13s_fix_recertify'
FIX = os.path.join(ROOT, 'D13s_work', 'fixed')
PY='<home>)/.venv/bin/python'
for cmd in [['step_18.py'], ['planner_sweep_ratio.py','1.5'], ['step_10_J192.py']]:
    r = subprocess.run([PY]+cmd, cwd=FIX, capture_output=True, text=True, timeout=3600)
    tail = r.stdout.strip().split(chr(10))
    print(f"--- {' '.join(cmd)} ---")
    for l in tail[-4:]: print("   ", l)
    if r.stderr: print("    STDERR:", r.stderr[-300:])
# and the big sweep JSON (already computed)
out = json.load(open(os.path.join(FIX,'sweep_cells_r1.0002.json')))
sup = max(o[2] for o in out if o[2] is not None)
print("--- sweep 1.0002 (from saved JSON): cells", len(out), " sup C =", sup)
