
import os, subprocess, time
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
PY='<home>)/.venv/bin/python'
t0 = time.time()
r = subprocess.run([PY, 'planner_sweep_ratio.py', '1.5'], cwd=FIX, capture_output=True, text=True, timeout=7200)
el = time.time()-t0
print(f"({el:.1f}s)")
print(r.stdout)
print(r.stderr[-1500:] if r.stderr else "")
with open(os.path.join(WORK,'progress.log'),'a') as f:
    f.write(f"D13s turn7: [100,400] ratio-1.05 rerun with fixed step_38: 29 cells, supC=194.1267 at (373.35,392.01) [was 194.3809 buggy].\n"
            f"D13s turn8: sweep ratio 1.5 rerun with fixed step_38 ({el:.0f}s):\n" + r.stdout.replace('\n','\n  ") + "\n")
