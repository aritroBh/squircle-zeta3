import subprocess
r = subprocess.run(["python3", "planner_check_naive.py"], cwd=base + "08_theory/COR5_corner_classicality/planner", capture_output=True, text=True)
print(r.stdout)
print(r.stderr[:500])
