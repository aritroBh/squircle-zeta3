import subprocess, os, time, datetime
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
with open(work+'/progress.log','a') as f:
    f.write(f'{datetime.datetime.now().isoformat()} LLL[w6x14] @3300 done (3200s): NO relation (spurious maxcoef ~2e19, res ~1e-3283). Next cumulative {14,15}.\n')
script=f'''#!/bin/bash
cd "{work}"
nohup nice -n 10 "{PY}" -B w6c_lll.py 3300 w6x14_15 > w6c_lll_x14_15_3300.out 2>&1 &
echo $! > w6c_lll_x14_15_3300.pid
'''
open(work+'/launch_x14_15.sh','w').write(script); os.chmod(work+'/launch_x14_15.sh',0o755)
subprocess.run(['bash',work+'/launch_x14_15.sh'])
time.sleep(45)
print('PID:',open(work+'/w6c_lll_x14_15_3300.pid').read())
print(subprocess.run(['bash','-lc',f'cat "{work}/w6c_lll_x14_15_3300.out"'],capture_output=True,text=True).stdout)
