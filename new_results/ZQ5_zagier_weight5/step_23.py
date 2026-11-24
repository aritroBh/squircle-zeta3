import sys, subprocess
src = r'''
import os
for p in ["../M4e_w5_lchi/planner/planner_check_Lrel.py",
          "../M4e_w5_lchi/planner/planner_eval_terms.py",
          "../M4i_w5_noL/verify_json.py"]:
    print("="*30, p, os.path.exists(p))
    if os.path.exists(p):
        t = open(p).read()
        print(t[:4000])
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',src],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-500:])
