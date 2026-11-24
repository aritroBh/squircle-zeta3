import os, subprocess, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
# remove the test pkl so it doesn't confuse
if os.path.exists('eval12_0.pkl'): os.remove('eval12_0.pkl')
if os.path.exists('eval12_0.log'): os.remove('eval12_0.log')
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_eval12.py','0','170'],
                      cwd=W, stdout=open(W+'/run_eval12_0.out','w'), stderr=subprocess.STDOUT)
open(W+'/run_eval12_0.pid','w').write(str(proc.pid))
print("launched chunk 0 (pts 0..169), pid", proc.pid, "at", time.strftime('%H:%M:%S'))
