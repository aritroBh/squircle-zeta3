import subprocess
PY='<project>/.venv/bin/python'
m5b='<project>/research/08_theory/M5b_w6_closed_form'
r=subprocess.run(['nice','-n','10',PY,'-B','w6_solve.py'],capture_output=True,text=True,cwd=m5b,timeout=1700)
print(r.stdout[-2000:]); print(r.stderr[-900:])
