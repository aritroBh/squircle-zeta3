import subprocess, os
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
cmd=f'cd {work} && nohup nice -n 10 {PY} -B m5b_precision.py w6 11500 11800 > w6_hp.out 2>&1 & echo $! > w6_hp.pid'
r=subprocess.run(['bash','-lc',cmd],capture_output=True,text=True)
print(r.stdout, r.stderr)
print(open(work+'/w6_hp.pid').read())
import time; time.sleep(5)
print(open(work+'/w6_hp.out').read())
