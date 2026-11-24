import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
with open('progress.log','a') as f:
    f.write(time.strftime('%H:%M:%S ')+"Step1d pipeline test (340 pts): m0 contributes EXACTLY 0 to A,B,Z,CLASS at all 340. Interpolated QA (deg8, 367 terms) from 255-pt rank-255 subset; D_A=0 on all 340. QA(1,1,1,1,1)=16472/25. Ansatz confirmed. Full 924-pt proof after campaign.\n")
print("logged")
time.sleep(400)
for st in [170,340,510,680,850]:
    lf=f'eval12_{st}.log'
    if os.path.exists(lf):
        print(lf, "->", open(lf).read().strip().split('\n')[-1])
print("time:", time.strftime('%H:%M:%S'))
