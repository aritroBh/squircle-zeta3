import subprocess
out = subprocess.run(["grep","-rln","dod_engine","../DOD1_dodecagon/"],capture_output=True,text=True).stdout
print(out)
out = subprocess.run(["grep","-rln","S_atom_decomp","../DOD1_dodecagon/","../DOD2_dodecagon_finish/"],capture_output=True,text=True).stdout
print(out)
