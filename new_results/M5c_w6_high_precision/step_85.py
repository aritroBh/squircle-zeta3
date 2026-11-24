import subprocess, os, time, datetime
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
with open(work+'/progress.log','a') as f:
    f.write(f'{datetime.datetime.now().isoformat()} LLL[w6x14_15] @3300 done (5770s): NO relation (spurious ~1e19, res ~1e-3283). Running full-lattice planted control @1300 before continuing ladder.\n')
script=f'''#!/bin/bash
cd "{work}"
nohup nice -n 10 "{PY}" -B w6c_control.py 1300 > w6c_control.out 2>&1 &
echo $! > w6c_control.pid
'''
open(work+'/launch_ctrl.sh','w').write(script); os.chmod(work+'/launch_ctrl.sh',0o755)
subprocess.run(['bash',work+'/launch_ctrl.sh'])
time.sleep(30)
print('PID:',open(work+'/w6c_control.pid').read())
print(subprocess.run(['bash','-lc',f'cat "{work}/w6c_control.out"'],capture_output=True,text=True).stdout)
