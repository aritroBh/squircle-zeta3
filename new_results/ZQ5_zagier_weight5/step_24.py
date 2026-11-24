import subprocess, sys
r = subprocess.run(['nice','-n','10',sys.executable,
                    '../M4e_w5_lchi/planner/planner_check_Lrel.py',
                    '../M4i_w5_noL/L_relation.json','830'],
                   capture_output=True,text=True,timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_11_theircheck: planner_check_Lrel at 830d\n")
