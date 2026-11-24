import os, subprocess, shutil
base='<project>/research/08_theory'
SRC=base+'/CHN4_p5_w11_locus'
W=base+'/CHN5_p5_w12_exact'
files=['point_engine.py','mzv_exact.py','kimilib.py','engine5.py','common.py','linalg_tools.py',
       'forms10.pkl','C4.txt','Q7p.txt','P7.txt','Q11.txt','Q11_poly.pkl',
       'R_a1a2a4.txt','R_factors.txt','subres1_A.txt','subres1_B.txt','gamma_final.log','gamma_final.pkl','gamma_final_pts.pkl']
for f in files:
    if os.path.exists(SRC+'/'+f):
        shutil.copy2(SRC+'/'+f, W+'/'+f)
    else:
        print("MISSING", f)
print(sorted(os.listdir(W)))
