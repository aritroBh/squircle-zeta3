import subprocess, os, time
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
script=f'''#!/bin/bash
cd "{work}"
nohup nice -n 10 "{PY}" -B w6_values.py 3500 > w6_vals_3500.out 2>&1 &
echo $! > w6_vals_3500.pid
'''
open(work+'/launch_vals.sh','w').write(script)
os.chmod(work+'/launch_vals.sh',0o755)
subprocess.run(['bash',work+'/launch_vals.sh'])
time.sleep(10)
print('PID:',open(work+'/w6_vals_3500.pid').read())
print(subprocess.run(['bash','-lc',f'tail -5 "{work}/w6_vals_3500.out"'],capture_output=True,text=True).stdout)
