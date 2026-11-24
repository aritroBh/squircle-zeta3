import os, subprocess, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_sign_nf.py'],
                      cwd=W, stdout=open(W+'/run_sign_nf.out','w'), stderr=subprocess.STDOUT)
proc.wait()
print(open(W+'/sign_nf.log').read())
print(open(W+'/run_sign_nf.out').read()[:2000])
