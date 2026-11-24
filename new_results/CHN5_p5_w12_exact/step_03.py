import subprocess, os
base='<project>/research/08_theory'
print(subprocess.run(['ls','-la',base+'/CHN4_p5_w11_locus'],capture_output=True,text=True).stdout)
print("=== CHN5 ===")
print(subprocess.run(['ls','-la',base+'/CHN5_p5_w12_exact'],capture_output=True,text=True).stdout)
