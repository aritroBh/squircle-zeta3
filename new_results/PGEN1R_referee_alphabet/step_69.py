import subprocess
out = subprocess.run(["grep","-rln","12-gon delta\|g2.*g3.*u\|data12\|DOD5","--include=step_*.py","../PGEN1_polygon_alphabet"],capture_output=True,text=True).stdout
print(out)
src = subprocess.run(["grep","-ln","12-gon delta","../PGEN1_polygon_alphabet/step_*.py"],capture_output=True,text=True).stdout
