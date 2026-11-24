
import os, shutil, subprocess, time
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
BUG = os.path.join(WORK, 'buggy')
os.makedirs(BUG, exist_ok=True)
PY='<home>)/.venv/bin/python'
shutil.copy(os.path.join(WORK,'D13g_all_N__step_38.py'), os.path.join(BUG,'step_38.py'))
shutil.copy(os.path.join(WORK,'D13h_final_cell__step_18.py'), os.path.join(BUG,'step_18.py'))
shutil.copy(os.path.join(WORK,'D13h_final_cell__planner_sweep_ratio.py'), os.path.join(BUG,'planner_sweep_ratio.py'))

# buggy [100,400]
t0=time.time()
r1 = subprocess.run([PY, 'step_18.py'], cwd=BUG, capture_output=True, text=True, timeout=1800)
last1 = [l for l in r1.stdout.strip().split(chr(10)) if 'sup C' in l]
print("BUGGY [100,400]:", last1, f"({time.time()-t0:.0f}s)")
# buggy sweep ratio 1.5
t0=time.time()
r2 = subprocess.run([PY, 'planner_sweep_ratio.py', '1.5'], cwd=BUG, capture_output=True, text=True, timeout=3600)
l2 = [l for l in r2.stdout.strip().split(chr(10)) if 'sup' in l or 'accepted' in l]
print("BUGGY sweep 1.5:", l2, f"({time.time()-t0:.0f}s)")
with open(os.path.join(WORK,'progress.log'),'a') as f:
    f.write(f"D13s turn11: buggy-code baselines: [100,400] {last1} ; sweep1.5 {l2}\n")
