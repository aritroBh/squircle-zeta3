import subprocess
PY='<project>/.venv/bin/python'
m5b='<project>/research/08_theory/M5b_w6_closed_form'
r=subprocess.run(['nice','-n','10',PY,'-B','g_eval5pp.py'],capture_output=True,text=True,cwd=m5b,timeout=1500)
print(r.stdout[-800:]); print(r.stderr[-600:])
