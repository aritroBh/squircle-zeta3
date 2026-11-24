import subprocess
r=subprocess.run(['nice','-n','10','<project>/.venv/bin/python','../M4e_w5_lchi/planner/planner_eval_R5nice.py','470'],capture_output=True,text=True,timeout=1200)
print(r.stdout)
print(r.stderr[-500:])
