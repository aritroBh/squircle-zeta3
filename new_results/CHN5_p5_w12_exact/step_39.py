import os, subprocess, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
print(subprocess.run(['bash','-c','ps aux | grep -i "run_eval12\|run_all" | grep -v grep'],capture_output=True,text=True).stdout)
print(subprocess.run(['ls','-la',W],capture_output=True,text=True).stdout[-2000:])
