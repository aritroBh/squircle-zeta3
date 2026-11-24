import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
time.sleep(200)
print(open('gamma_trace.log').read()[-300:])
print("OUT:", open('run_gamma_trace.out').read()[:300])
print("gamma_trace.pkl exists:", os.path.exists('gamma_trace.pkl'))
