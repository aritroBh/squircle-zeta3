import subprocess, os
base='<project>/research/08_theory'
SRC=base+'/CHN4_p5_w11_locus'
# search step files for det11 construction and rank verification
out=subprocess.run(['grep','-l','orbits11\|det11\|rank'], [SRC+'/'+f for f in os.listdir(SRC) if f.startswith('step_')],capture_output=True,text=True).stdout
print(out)
