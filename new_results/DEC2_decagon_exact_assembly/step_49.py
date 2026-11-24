import subprocess
out = subprocess.run(['grep','-rl','Jtab', '<project>/research/08_theory/K1_nonlocal_law'], capture_output=True, text=True).stdout
print(out)
base = "<project>/research/08_theory/K1_nonlocal_law"
import os
if os.path.isdir(base):
    for f in sorted(os.listdir(base)):
        if f.endswith('.py') and 'Jtab' in open(os.path.join(base,f)).read():
            txt = open(os.path.join(base,f)).read()
            if 'pickle.dump' in txt:
                print("=====", f)
                print(txt[:3000])
                break
