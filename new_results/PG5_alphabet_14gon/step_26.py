import subprocess
out = subprocess.run(["grep","-rln","E3","../DOD1_dodecagon/","--include=step_*.py"],capture_output=True,text=True).stdout
print(out)
out = subprocess.run(["grep","-ln","S_atom_decomp\|atom_coefs","../DOD1_dodecagon/step_*.py","../DOD2_dodecagon_finish/*.py"],capture_output=True,text=True).stdout
print(out)
