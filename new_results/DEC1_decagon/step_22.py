import subprocess
d='research/08_theory/OCT6_octagon_closed/'
out=subprocess.run(['grep','-ln','MixedW\|mixed',*[d+f'step_{i:02d}.py' for i in range(6,64)]],capture_output=True,text=True)
print(out.stdout)
out=subprocess.run(['sed','-n','1,80p',d+'progress.log'],capture_output=True,text=True)
print(out.stdout)