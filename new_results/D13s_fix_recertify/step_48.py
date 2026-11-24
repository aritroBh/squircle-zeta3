
import json, os, subprocess, time
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
BUG = os.path.join(WORK, 'buggy')
PY='<home>)/.venv/bin/python'
out = json.load(open(os.path.join(BUG, 'sweep_cells_r1.5.json')))
sup = max(o[2] for o in out if o[2] is not None)
argm = [o[:2] for o in out if o[2]==sup][0]
print(f"BUGGY ratio-1.5 sweep: sup C = {sup:.6f} at {argm}")

t0=time.time()
r = subprocess.run([PY, 'planner_sweep_ratio.py', '1.0002'], cwd=BUG, capture_output=True, text=True, timeout=3600)
el = time.time()-t0
print(f"BUGGY 1.0002 sweep ({el:.0f}s):")
print(r.stdout)
print(r.stderr[-800:] if r.stderr else "")
out2 = json.load(open(os.path.join(BUG, 'sweep_cells_r1.0002.json')))
sup2 = max(o[2] for o in out2 if o[2] is not None)
argm2 = [o[:2] for o in out2 if o[2]==sup2][0]
print(f"BUGGY ratio-1.0002 sweep: sup C = {sup2:.6f} at {argm2}, cells {len(out2)}")
with open(os.path.join(WORK,'progress.log'),'a') as f:
    f.write(f"D13s turn12: buggy sweep1.5 exact sup={sup:.6f} at {argm}; buggy sweep1.0002 sup={sup2:.6f} at {argm2} ({len(out2)} cells, {el:.0f}s).\n")
