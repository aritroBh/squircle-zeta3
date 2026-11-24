import os, shutil, sys, time
base='<project>/research/08_theory'
SRC=base+'/CHN4_p5_w11_locus'; W=base+'/CHN5_p5_w12_exact'
shutil.copy2(SRC+'/survey_lib.py', W+'/survey_lib.py')
os.chdir(W); sys.path.insert(0,W)
import numpy as np, importlib
import survey_lib as SL
importlib.reload(SL)
# test V10 points at a1=0.01, a2=2.74 (near the known Gamma point)
pts=SL.v10_points(0.01, 2.74)
print("V10 points at a1=0.01,a2=2.74:", pts)
# Q11 (=-P7/10) sign there:
for (a4v,a5v) in pts:
    print(f"  a4={a4v:.6f} a5={a5v:.6f}  Q11={SL.q11_np([0.01,2.74,1,a4v,a5v]):+.6g}")
