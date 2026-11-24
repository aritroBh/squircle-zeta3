import subprocess, time
PY='<project>/.venv/bin/python'
m5b='<project>/research/08_theory/M5b_w6_closed_form'
t0=time.time()
r=subprocess.run(['nice','-n','10',PY,'-B','m5b_precision.py','w6','5800'],capture_output=True,text=True,cwd=m5b,timeout=2700)
print('w6@5800 elapsed', time.time()-t0); print(r.stdout[-600:]); print(r.stderr[-500:])
