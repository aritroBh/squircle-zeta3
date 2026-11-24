import os, subprocess, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_factor.py'],
                      cwd=W, stdout=open(W+'/run_factor.out','w'), stderr=subprocess.STDOUT)
open(W+'/run_factor.pid','w').write(str(proc.pid))
print("launched factor pid", proc.pid, "at", time.strftime('%H:%M:%S'))
time.sleep(150)
print(open('factor.log').read() if os.path.exists('factor.log') else "(no log yet)")
