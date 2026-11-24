import subprocess, time, os
PY='<project>/.venv/bin/python'
m5b='<project>/research/08_theory/M5b_w6_closed_form'
t0=time.time()
r=subprocess.run(['nice','-n','10',PY,'-B','m5b_precision.py','w6','8000'],capture_output=True,text=True,cwd=m5b,timeout=1700)
print('w6@8000 elapsed', time.time()-t0)
out=r.stdout
i=out.find('"digits"')
print(out[i-60:i+250] if i>0 else out[-300:]); print(r.stderr[-200:])
