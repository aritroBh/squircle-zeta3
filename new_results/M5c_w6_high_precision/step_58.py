import subprocess, os, time, datetime
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
with open(work+'/progress.log','a') as f:
    f.write(f'{datetime.datetime.now().isoformat()} LLL @3000 done (3250s, load~280): NO target relation. Best spurious rows maxcoef 4.9e17-8.2e17, residuals ~1e-2984. Genuine coeffs must exceed ~1e17.3 -> run 3300.\n')
script=f'''#!/bin/bash
cd "{work}"
nohup nice -n 10 "{PY}" -B w6c_lll.py 3300 > w6c_lll_3300.out 2>&1 &
echo $! > w6c_lll_3300.pid
'''
open(work+'/launch_lll3300.sh','w').write(script); os.chmod(work+'/launch_lll3300.sh',0o755)
subprocess.run(['bash',work+'/launch_lll3300.sh'])
time.sleep(45)
print('PID:',open(work+'/w6c_lll_3300.pid').read())
print(subprocess.run(['bash','-lc',f'cat "{work}/w6c_lll_3300.out"'],capture_output=True,text=True).stdout)
