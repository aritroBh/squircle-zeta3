import subprocess, os
wt3 = "<project>/research/08_theory/WT3_weight3_all_k/"
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","verify_wt3.py"],
                 cwd=wt3,capture_output=True,text=True,env=env,timeout=1700)
print(r.stdout)
print("STDERR (tail):", r.stderr[-1500:])
work = "<project>/research/08_theory/TOW1R_review/"
with open(work+"progress.log","a") as fh: fh.write("step_24: reran authors' verify_wt3.py fully\n")
