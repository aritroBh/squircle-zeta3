import subprocess, os, time
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
script=f'''#!/bin/bash
cd "{work}"
nohup nice -n 10 "{PY}" -B w6c_lll.py 3300 w6x14_15_9 > w6c_lll_x14_15_9_3300.out 2>&1 &
echo $! > w6c_lll_x14_15_9_3300.pid
'''
open(work+'/launch_x14_15_9.sh','w').write(script); os.chmod(work+'/launch_x14_15_9.sh',0o755)
subprocess.run(['bash',work+'/launch_x14_15_9.sh'])
time.sleep(45)
print('PID:',open(work+'/w6c_lll_x14_15_9_3300.pid').read())
print(subprocess.run(['bash','-lc',f'cat "{work}/w6c_lll_x14_15_9_3300.out"'],capture_output=True,text=True).stdout)
