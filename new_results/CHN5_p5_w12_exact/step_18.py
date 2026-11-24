import subprocess, os
base='<project>/research/08_theory'
SRC=base+'/CHN4_p5_w11_locus'
out=subprocess.run(['grep','-l','orbits11']+[os.path.join(SRC,f) for f in os.listdir(SRC) if f.startswith('step_')],capture_output=True,text=True).stdout
print("orbits11 files:", out)
out=subprocess.run(['grep','-l','rank']+[os.path.join(SRC,f) for f in os.listdir(SRC) if f.startswith('step_')],capture_output=True,text=True).stdout
print("rank files:", out)
