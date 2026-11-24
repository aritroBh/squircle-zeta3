import time, os
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
time.sleep(300)
print(open('eval12_0.log').read()[-200:])
