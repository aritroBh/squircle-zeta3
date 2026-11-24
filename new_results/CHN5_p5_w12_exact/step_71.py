import os, subprocess, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
if os.path.exists('gamma_trace.log'): os.remove('gamma_trace.log')
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_gamma_trace.py'],
                      cwd=W, stdout=open(W+'/run_gamma_trace.out','w'), stderr=subprocess.STDOUT)
open(W+'/run_gamma_trace.pid','w').write(str(proc.pid))
print("launched gamma trace pid", proc.pid, "at", time.strftime('%H:%M:%S'))
time.sleep(120)
print(open('gamma_trace.log').read() if os.path.exists('gamma_trace.log') else "(no log yet)")
print(open('run_gamma_trace.out').read()[:500])
