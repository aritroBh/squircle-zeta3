import subprocess, os
tow = "<project>/research/08_theory/TOW1_tower_law/"
env = dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r = subprocess.run(["nice","-n","19","<project>/.venv/bin/python","verify_tower.py"],
                   cwd=tow, capture_output=True, text=True, env=env, timeout=1700)
print(r.stdout)
print("STDERR (tail):", r.stderr[-2000:])
work = "<project>/research/08_theory/TOW1R_review/"
with open(work+"progress.log","a") as fh: fh.write("step_04: reran authors' verify_tower.py in place\n")
