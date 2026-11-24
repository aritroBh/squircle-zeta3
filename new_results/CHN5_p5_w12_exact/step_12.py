import os, subprocess
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
for f in ['run.log','step_01.py','step_02.py','step_03.py']:
    p=W+'/'+f
    if os.path.exists(p):
        print("="*20,f,"="*20)
        print(open(p).read())
print("==== CHN5 brief ====")
print(open(base+'/CHN5_p5_w12_exact.brief').read()[:4000] if os.path.exists(base+'/CHN5_p5_w12_exact.brief') else "none")
