import sys, subprocess
print("python tool executable:", sys.executable)
VENV='<project>/.venv/bin/python'
r = subprocess.run([VENV,'-c','import sys,mpmath,sympy,flint;print("venv ok",sys.version.split()[0]);import mpmath as mp;print("mpmath",mp.__version__)'],capture_output=True,text=True)
print(r.stdout, r.stderr)
# set up workdir
import os
WD='<project>/research/08_theory/FLD1_field_rule_proof'
os.chdir(WD)
# copy the template modules we need
import shutil
B='<project>/research/10_frontier/W1_corner_family/B_taylor_pi4/'
for f in ['zsum.py','hpl.py']:
    shutil.copy(B+f, WD+'/'+f)
print("workdir contents:", os.listdir(WD))
with open(WD+'/progress.log','a') as fh:
    fh.write("setup: copied zsum.py hpl.py; environment verified\n")
print("done")