import subprocess, os, time
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
script=f'''#!/bin/bash
cd "{work}"
nohup nice -n 10 "{PY}" -B m5b_precision.py w6 11500 11800 > w6_hp.out 2>&1 &
echo $! > w6_hp.pid
'''
open(work+'/launch_w6.sh','w').write(script)
os.chmod(work+'/launch_w6.sh',0o755)
r=subprocess.run(['bash',work+'/launch_w6.sh'],capture_output=True,text=True)
print(r.stdout, r.stderr)
time.sleep(8)
print('PID:',open(work+'/w6_hp.pid').read())
print(open(work+'/w6_hp.out').read())
