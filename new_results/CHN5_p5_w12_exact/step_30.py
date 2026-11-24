import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
print(open('eval12_0.log').read() if os.path.exists('eval12_0.log') else "no log yet")
print("time now:", time.strftime('%H:%M:%S'))
