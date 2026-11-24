import os, subprocess
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
print(subprocess.run(['bash','-c',f'ls -la {W} | grep -E "eval12|run_all|run_eval12"'],capture_output=True,text=True).stdout)
print("run_all.out:", repr(open(W+'/run_all.out').read()) if os.path.exists(W+'/run_all.out') else "missing")
