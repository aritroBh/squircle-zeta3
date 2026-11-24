import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
time.sleep(580)
for st in [170,340,510,680,850]:
    lf=f'eval12_{st}.log'
    if os.path.exists(lf):
        print(lf, "->", open(lf).read().strip().split('\n')[-1])
print("time:", time.strftime('%H:%M:%S'))
