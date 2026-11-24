import subprocess
PY='<project>/.venv/bin/python'
A2='<project>/research/08_theory/ALPHA2_angle_leading_law'
out = open(f'{A2}/kern18_new.out','w')
proc = subprocess.Popen(['nice','-n','10',PY,f'{A2}/kern18.py','new'], stdout=out, stderr=subprocess.STDOUT)
print("launched new", proc.pid)
with open(f'{A2}/progress.log','a') as fh:
    fh.write('step: F_d extracted d<=7 (8 angles); leading law verified 8 angles x w=2..7 (44-51 digits); Master Lemma corner version confirmed; kern18 old done, new launched\n')
