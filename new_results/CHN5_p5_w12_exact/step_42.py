import os, subprocess, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
time.sleep(40)
print(subprocess.run(['bash','-c','ps aux | grep run_eval12 | grep -v grep | head -3'],capture_output=True,text=True).stdout)
print("eval12_170.log:", open(W+'/eval12_170.log').read() if os.path.exists(W+'/eval12_170.log') else "(not yet)")
