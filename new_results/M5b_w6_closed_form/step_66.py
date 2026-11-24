import subprocess
PY='<project>/.venv/bin/python'
m5b='<project>/research/08_theory/M5b_w6_closed_form'
r=subprocess.run(['nice','-n','10',PY,'-B','w5pp_lll.py','1300'],capture_output=True,text=True,cwd=m5b,timeout=1700)
out=r.stdout
# print only the summary parts
import re
print(out[:600])
i=out.find('"rels"')
print(out[i:i+400] if i>0 else out[-1500:])
print(r.stderr[-500:])
