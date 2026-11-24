import subprocess, os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'-B','m5b_precision.py','w6','600'],capture_output=True,text=True,cwd=m5b,timeout=900)
print(r.stdout[-2000:]); print(r.stderr[-1000:])
