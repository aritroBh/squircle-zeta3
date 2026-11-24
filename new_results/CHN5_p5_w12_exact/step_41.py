import os, pickle
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
print("eval12_0.pkl exists:", os.path.exists(W+'/eval12_0.pkl'), os.path.getsize(W+'/eval12_0.pkl') if os.path.exists(W+'/eval12_0.pkl') else 0)
runner = '''#!/bin/bash
cd "<project>/research/08_theory/CHN5_p5_w12_exact"
PY=<project>/.venv/bin/python
for start in 170 340 510 680 850; do
  cnt=170
  if [ $start -eq 850 ]; then cnt=74; fi
  nice -n 10 "<project>/.venv/bin/python" run_eval12.py $start $cnt >> run_all.log 2>&1
done
echo "ALL CHUNKS DONE" >> run_all.log
'''
open(W+'/run_all.sh','w').write(runner)
import subprocess, time
proc=subprocess.Popen(['nice','-n','10','bash','run_all.sh'],
                      cwd=W, stdout=open(W+'/run_all.out','w'), stderr=subprocess.STDOUT)
open(W+'/run_all.pid','w').write(str(proc.pid))
print("relaunched run_all.sh pid", proc.pid, "at", time.strftime('%H:%M:%S'))
