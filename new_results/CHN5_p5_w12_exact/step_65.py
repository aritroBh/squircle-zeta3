import os, subprocess, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_interp12.py'],
                      cwd=W, stdout=open(W+'/run_interp12.out','w'), stderr=subprocess.STDOUT)
open(W+'/run_interp12.pid','w').write(str(proc.pid))
print("launched interp pid", proc.pid, "at", time.strftime('%H:%M:%S'))
time.sleep(90)
print(open(W+'/interp12.log').read() if os.path.exists(W+'/interp12.log') else "(no log yet)")
print(open(W+'/run_interp12.out').read())
