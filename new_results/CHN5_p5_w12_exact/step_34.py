import os, subprocess, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
runner = r'''#!/bin/bash
cd <project>/research/08_theory/CHN5_p5_w12_exact
PY=<project>/.venv/bin/python
for start in 170 340 510 680 850; do
  cnt=170
  if [ $start -eq 850 ]; then cnt=74; fi
  nice -n 10 $PY run_eval12.py $start $cnt >> run_all.log 2>&1
done
echo "ALL CHUNKS DONE" >> run_all.log
'''
open(W+'/run_all.sh','w').write(runner)
os.chmod(W+'/run_all.sh',0o755)
proc=subprocess.Popen(['nice','-n','10','bash',W+'/run_all.sh'],
                      cwd=W, stdout=open(W+'/run_all.out','w'), stderr=subprocess.STDOUT)
open(W+'/run_all.pid','w').write(str(proc.pid))
print("launched run_all.sh pid", proc.pid, "at", time.strftime('%H:%M:%S'))
# quick sanity of chunk-0 data: load eval12_0.pkl and inspect one point's m0
import pickle
d0=pickle.load(open('eval12_0.pkl','rb'))
print("chunk0 pts:", len(d0))
from fractions import Fraction as F
ex=d0[0]
print("pt:", ex['pt'])
print("m0 coords 7..11:", ex['m0'][7:])
print("m0 CLASS:", ex['m0'][9]+F(14,9)*ex['m0'][10])
