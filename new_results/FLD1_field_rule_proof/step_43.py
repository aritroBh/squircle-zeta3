import subprocess, time
WD='<project>/research/08_theory/FLD1_field_rule_proof'
VENV='<project>/.venv/bin/python'
t=time.time()
r = subprocess.run(['nice',VENV,'step_25_construct.py','5'],capture_output=True,text=True,cwd=WD)
print("time", round(time.time()-t,1),"s")
for line in r.stdout.splitlines():
    if 'residual' in line or 'MISMATCH' in line: print(line[:250]); continue
    print(line)
print(r.stderr[-3000:] if r.returncode else '')