import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
time.sleep(200)
print(open('eval12_850.log').read().strip().split('\n')[-1])
print("run_all.log tail:", open('run_all.log').read()[-60:])
