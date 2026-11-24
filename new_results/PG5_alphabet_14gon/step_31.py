import subprocess
out = subprocess.run(["grep","-ln","S_atom_decomp\|dod_engine\|import dod_engine\|from dod_engine","../DOD1_dodecagon/step_*.py"],capture_output=True,text=True).stdout
print(out)
