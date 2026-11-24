import subprocess, time
PY='<project>/.venv/bin/python'
A2='<project>/research/08_theory/ALPHA2_angle_leading_law'
out = open(f'{A2}/verify_alpha2.out','w')
proc = subprocess.Popen(['nice','-n','10',PY,f'{A2}/verify_alpha2.py'], stdout=out, stderr=subprocess.STDOUT)
print("launched", proc.pid)
time.sleep(90)
print(open(f'{A2}/verify_alpha2.out').read())
