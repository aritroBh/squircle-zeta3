import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
time.sleep(540)
print(open('gamma_trace.log').read()[-2200:])
