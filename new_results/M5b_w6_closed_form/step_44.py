import os, json, glob
m5b='<project>/research/08_theory/M5b_w6_closed_form'
for f in sorted(glob.glob(m5b+'/w5pp*')):
    print(os.path.basename(f), os.path.getsize(f))
print(open(m5b+'/progress.log').read())
import subprocess
print(subprocess.run(['pgrep','-fl','m5b_precision'],capture_output=True,text=True).stdout)
