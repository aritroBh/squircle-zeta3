import shutil, subprocess, os
m4c='<project>/research/08_theory/M4c_w5_closed_form'
m5b='<project>/research/08_theory/M5b_w6_closed_form'
shutil.copy(f'{m4c}/w5_precision.py', f'{m5b}/w5_precision.py')
print('copied w5_precision.py')
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'-B','w5pp_cobracket.py'],capture_output=True,text=True,cwd=m5b,timeout=1700)
print(r.stdout[-3000:]); print(r.stderr[-900:])
