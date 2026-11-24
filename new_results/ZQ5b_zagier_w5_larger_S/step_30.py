import subprocess
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S') if False else None
import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
print(subprocess.run(['grep','-n','unit_emb','w5blib.py'],capture_output=True,text=True,cwd='.').stdout)
print(subprocess.run(['sed','-n','108,118p','w5blib.py'],capture_output=True,text=True,cwd='.').stdout)
