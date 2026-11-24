import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
time.sleep(200)
print(open('gamma_trace2.log').read()[-400:])
print("OUT:", open('run_gamma_trace2.out').read()[:200])
