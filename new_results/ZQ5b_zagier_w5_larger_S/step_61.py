import os, subprocess, time
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
# launch in background
os.system('rm -f step09b.log')
proc = subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','-u','step_09b_d13_kernel.py'],
                        stdout=open('step09b.log','w'), stderr=subprocess.STDOUT, cwd='.')
print("launched pid", proc.pid)
time.sleep(60)
print(open('step09b.log').read())
print("--- still running:", proc.poll() is None)
