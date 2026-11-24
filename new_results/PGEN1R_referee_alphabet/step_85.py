import subprocess
out = subprocess.run(["grep","-rln","L3(xi)\|xi_e\|single-valued\|L_3","--include=step_*.py","../PGEN1_polygon_alphabet"],capture_output=True,text=True).stdout
print(out[:600])
src = subprocess.run(["sed","-n","1,60p","../PGEN1_polygon_alphabet/step_82.py"],capture_output=True,text=True).stdout
print(src)
