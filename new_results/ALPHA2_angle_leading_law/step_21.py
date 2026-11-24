import subprocess, os
PY='<project>/.venv/bin/python'
A2='<project>/research/08_theory/ALPHA2_angle_leading_law'
out = open(f'{A2}/kern_new.out','w')
proc = subprocess.Popen(['nice','-n','10',PY,f'{A2}/kern_new.py'], stdout=out, stderr=subprocess.STDOUT)
print("launched pid", proc.pid)
open(f'{A2}/kern_new.pid','w').write(str(proc.pid))
