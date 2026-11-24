import subprocess, os, time, datetime
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
script=f'''#!/bin/bash
cd "{work}"
nohup nice -n 10 "{PY}" -B w6c_lll.py 3000 > w6c_lll_3000.out 2>&1 &
echo $! > w6c_lll_3000.pid
'''
open(work+'/launch_lll3000.sh','w').write(script); os.chmod(work+'/launch_lll3000.sh',0o755)
subprocess.run(['bash',work+'/launch_lll3000.sh'])
time.sleep(60)
print('PID:',open(work+'/w6c_lll_3000.pid').read())
print(subprocess.run(['bash','-lc',f'cat "{work}/w6c_lll_3000.out"'],capture_output=True,text=True).stdout)
