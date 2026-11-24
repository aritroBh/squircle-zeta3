import subprocess, time
WD='<project>/research/08_theory/FLD1_field_rule_proof'
VENV='<project>/.venv/bin/python'
# fast numeric mu table k=6 (N=40)
t=time.time()
r = subprocess.run(['nice',VENV,'step_21_fastmu.py','6','40'],capture_output=True,text=True,cwd=WD)
print("k=6 fast table time", round(time.time()-t,1),"s")
for line in r.stdout.splitlines():
    if 'mu_' in line or 'computed' in line or 'matches' in line or 'saved' in line: print(line)
print(r.stderr[-500:] if r.returncode else '')